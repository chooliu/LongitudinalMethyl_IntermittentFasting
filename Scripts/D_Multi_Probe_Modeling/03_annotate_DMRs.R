# ==============================================================================
# D/03_annotate_DMRs.R
# summarize gene overlap, DMPs in each candidate DMR
# ==============================================================================



# usually run interactively, ~3 hours ------------------------------------------

screen -S combpprocessing
qlogin -R rusage[mem=16]

module load R/4.1.0_beta
module load gcc/7.4.0 
cd $HOME/Analyses/DRIFT2_Methylation/
R



# load gene features, sesame annots --------------------------------------------


library(tidyverse)
library(data.table)
library(annotatr)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(sesameData)
library(writexl)

EPIC.hg38.manifest <-
  sesameDataGet("EPIC.hg38.manifest")

genewise_annotations <-
  build_annotations(
    genome = "hg38",
    annotations = "hg38_basicgenes"
  )

load("Output/20220425/probe_annotations.Rdata")




# annotate each candidate DMR --------------------------------------------------
# assumes that maintains similar filename as DMPs from part C

summarize_DMR_results <-
  function(dmp_path, dmr_path, padding = 0) {
    
    # read DMP and DMR results files
    dmp_results <-
      fread(dmp_path) %>%
      left_join(., probe_annotations, by = "Probe")
    dmr_results <- fread(dmr_path)
    
    if (nrow(dmr_results) == 0) {
      return(NULL)
    }
    
    # subset genome annotations to DMR overlap,
    # then compile gene names to format "GENEX (3' UTR/exon)"
    dmr_anno <- 
      annotate_regions(
        regions = GRanges(seqnames = dmr_results$`#chrom`,
                          ranges = IRanges(
                            start = (dmr_results$start - padding),
                            end = dmr_results$end + padding)),
        annotations = genewise_annotations
      ) %>%
      as.data.frame(.) %>%
      mutate(start = start + padding,
             end = end - padding)
    
    dmr_anno_summary <-
      dmr_anno %>%
      mutate(group =
               paste0(seqnames, start, end, annot.symbol, annot.type),
             symbol_fxn = paste0(annot.symbol, " (",
                                 ( gsub("hg38_genes|_", "", annot.type) %>%
                                     gsub("s", "", .) ),
                                 ")")) %>%
      filter(!duplicated(group)) %>%
      group_by(seqnames, start, end, annot.symbol) %>%
      summarize(types = unique(annot.type) %>% gsub("hg38_genes|_", "", .) %>%
                  gsub("s", "", .) %>% paste0(collapse = "/")) %>%
      ungroup() %>%
      mutate(symbol_types = paste0(annot.symbol, " (", types, ")")) %>%
      group_by(seqnames, start, end) %>% 
      summarize(genes = unique(symbol_types) %>% paste(collapse = "; ") %>%
                  gsub("NA \\(", "N/A \\(", .)) %>%
      set_names(c("DMR Chromosome", "DMR Start", "DMR End", "DMR Annotation"))
    
    # intersect DMR with DMP list
    dmr_to_dmp_mapping <-
      dmr_results %>%
      apply(., 1,
            function(row) {
              probes_in_region <-
                dmp_results[(Chr == gsub("chr", "", row[1])) &
                              between(Start, row[2], row[3]), ] %>%
                arrange(bacon_Pz)
              out <-
                probes_in_region %>%
                dplyr::slice(1) %>%
                transmute(`Top DMP` = Probe,
                          `DMP Position` = Start,
                          `DMP Effect` = Estimate,
                          `DMP StdErr` = SE,
                          `DMP P-Value` = bacon_Pz,
                          `DMP FDR` = fdr,
                          `DMP Gene` = Gene,
                          BetaConversion1,
                          BetaConversion2)
              # do DMPs in the DMR have same direction of effect
              bind_cols(out,
                        `# Probes Within Region` = nrow(probes_in_region),
                        `Probe Directions (+/-)` = paste0(sum(probes_in_region$Estimate > 0), "/", sum(probes_in_region$Estimate < 0)),
                        `Proportion +` = sum(probes_in_region$Estimate > 0) / nrow(probes_in_region))
            })
    
    # join combp & DMP output, export
    dmr_results_final <-
      dmr_results %>%
      transmute(`DMR Chromosome` = `#chrom`,
                `DMR Start` = start,
                `DMR End` = end,
                `DMR P-value` = z_p,
                `DMR Adjusted P-Value` = z_sidak_p,
                `DMR Signif` = if_else(`DMR Adjusted P-Value` < 0.10, "*", "")) %>%
      left_join(dmr_anno_summary, by = c("DMR Chromosome", "DMR Start", "DMR End")) %>%
      mutate(`DMR Annotation` = if_else(is.na(`DMR Annotation`), "", `DMR Annotation`)) %>%
      bind_cols(dmr_to_dmp_mapping %>% bind_rows())  %>%
      arrange(`DMR P-value`) %>%
      mutate(`DMR Chromosome` = gsub("chr", "", `DMR Chromosome`)) %>%
      mutate_at(.vars = c("DMR Start", "DMR End", "DMP Position"), as.character) %>%
      filter(`DMR Signif` == "*") # for this manuscript, export only signif DMRs
    
    dmr_results_final
  }

# excel worksheet name
convert_filepath_to_outname <-
  function(filepath) {
    basename(filepath) %>%
      gsub(".regions-t.bed", "", .) %>%
      gsub("Methyl0", "M0", .) %>%
      gsub("cellcomp", "cellc", .) %>%
      gsub("covar", "set", .) %>%
      gsub("RUV", "", .) %>%
      gsub("__", "_", .)
  }

# same z-value <--> beta conversions as DMP annotation
convert_effect_sizes <-
  function(df, model_type) {
    if (is.null(df)) {
      return(NULL)
    }
    if (model_type == 1) {
      df_out <-
        df %>%
        mutate(`Beta Conversion` = BetaConversion1,
               `DMP Effect Beta` = `DMP Effect` / BetaConversion1) %>%
        dplyr::select(-BetaConversion1, -BetaConversion2)
    } else {
      df_out <-
        df %>%
        mutate(`Beta Conversion` = BetaConversion2,
               `DMP Effect Beta` = `DMP Effect` / BetaConversion2) %>%
        dplyr::select(-BetaConversion1, -BetaConversion2)
    }
    return(df_out)
  }



# run annotations and export (baseline only models) ----------------------------

dmr_path <-
  list.files("./Output/20220425/combp/results", full.names = T) %>%
  .[grepl("regions-t.bed", .)] %>%
  .[grepl("BLonly", .)]
dmp_path <-
  gsub("combp/results", "ResultsTables", dmr_path) %>%
  gsub(".regions-t.bed", ".tsv", .)

compiled_comp_results <-
  tibble(dmr_path = dmr_path,
       dmp_path = dmp_path) %>%
  pmap(summarize_DMR_results)

compiled_comp_results <-
  map(1:length(compiled_comp_results),
      function(i) {
        convert_effect_sizes(
          df = compiled_comp_results[[i]],
          model_type = if_else(grepl("RUVcellc", dmr_path[i]), 2, 1)) })

excel_names <-
  dmr_path %>% sapply(convert_filepath_to_outname)
write_xlsx(
  compiled_comp_results %>%
    set_names(gsub("BLonly_", "", excel_names)) %>% 
    .[map_lgl(., ~ !(is.null(.x) | nrow(.x) == 0))] %>%
    map(as_tibble),
  "Output/20220425/CompiledResults/DMRs_BLonly.xlsx",
  format_headers = F
)



# baseline methylation --> change models ---------------------------------------

dmr_path <-
  list.files("./Output/20220425/combp/results", full.names = T) %>%
  .[grepl("regions-t.bed", .)] %>%
  .[grepl("Methyl0_", .)] %>% .[grepl("_30_|_60_", .)]
dmp_path <-
  gsub("combp/results", "ResultsTables", dmr_path) %>%
  gsub(".regions-t.bed", ".tsv", .)

compiled_comp_results <-
  tibble(dmr_path = dmr_path,
         dmp_path = dmp_path) %>%
  pmap(summarize_DMR_results)

compiled_comp_results <-
  map(1:length(compiled_comp_results),
      function(i) {
        convert_effect_sizes(
          df = compiled_comp_results[[i]],
          model_type = if_else(grepl("RUVcellc", dmr_path[i]), 2, 1)) })

excel_names <-
  dmr_path %>% sapply(convert_filepath_to_outname)
write_xlsx(
  compiled_comp_results %>%
    set_names(excel_names) %>% 
    .[map_lgl(., ~ is.null(.x) %>% `!`)] %>%
    map(as_tibble),
  "Output/20220425/CompiledResults/DMRs_ShortTermChange.xlsx",
  format_headers = F
)




# time/treatment --> methylation models ----------------------------------------

dmr_path <-
  list.files("./Output/20220425/combp/results", full.names = T) %>%
  .[grepl("regions-t.bed", .)] %>%
  .[grepl("Trt|Time", .)]
dmp_path <-
  gsub("combp/results", "ResultsTables", dmr_path) %>%
  gsub(".regions-t.bed", ".tsv", .)

compiled_comp_results <-
  tibble(dmr_path = dmr_path,
         dmp_path = dmp_path) %>%
  pmap(summarize_DMR_results)

compiled_comp_results <-
  map(1:length(compiled_comp_results),
      function(i) {
        convert_effect_sizes(
          df = compiled_comp_results[[i]],
          model_type = if_else(grepl("RUVcellc", dmr_path[i]), 2, 1)) }) %>%
  map(~ mutate(.x, EffectBeta = Effect * BetaConversion) )

excel_names <-
  dmr_path %>% sapply(convert_filepath_to_outname)
write_xlsx(
  compiled_comp_results %>%
    set_names(excel_names) %>% 
    .[map_lgl(., ~ is.null(.x) %>% `!`)] %>%
    map(as_tibble),
  "Output/20220425/CompiledResults/DMRs_TimeTreatment.xlsx",
  format_headers = F
)


