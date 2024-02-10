# ==============================================================================
# D/07_gometh.R
# pathway enrichment via gometh
# ==============================================================================



screen -S gomethTime
qlogin -R rusage[mem=16]

module load R/4.1.0_beta
cd $HOME/Analyses/DRIFT2_Methylation/
R


library(data.table)
library(tidyverse)
library(readxl)
library(GenomicRanges)
library(writexl)
library(missMethyl)

dir.create("Output/20220425/GOMeth/")



ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)



convert_filepath_to_outname <-
  function(filepath) {
    basename(filepath) %>%
      gsub(".tsv", "", .) %>%
      gsub("Methyl0", "M0", .) %>%
      gsub("cellcomp", "cellc", .) %>%
      gsub("covar", "set", .) %>%
      gsub("RUV", "", .) 
  }


run_enrichment_test <-
  function(filepath_in, genomicfeaturestest) {
    
    results_table <-
      fread(filepath_in) %>%
      filter(!is.na(bacon_Pz))
    
    num_probes <-
      max(500, results_table$bacon_Pz < 1e-4) 
    
    gometh_results <-
      gometh(sig.cpg = results_table$Probe[1:num_probes],
             all.cpg = results_table$Probe,
             collection	= 'KEGG',
             array.type = "EPIC",
             prior.prob = T,
             genomic.features = genomicfeaturestest,
             anno = ann,
             sig.genes = T) %>%
      as_tibble(rownames = "KEGGID")
    
    gometh_results <-
      gometh_results %>% filter(N >= 10) %>% arrange(P.DE)
    gometh_results$FDR <- p.adjust(gometh_results$P.DE, "fdr")
    
    return(gometh_results)
    
  }





# short term change ------------------------------------------------------------

filepaths <-
  list.files(path = "Output/20220425/ResultsTables", full.names = T) %>%
  .[grepl("_30_|_60_", .)] %>%
  .[grepl("covarB", .)]
excel_names <-
  sapply(filepaths, convert_filepath_to_outname)

enrichment_results_promoter <-
  map(filepaths, run_enrichment_test,
      genomicfeaturestest = c("TSS1500", "TSS200", "1stExon"))

write_xlsx(enrichment_results_promoter %>% 
             set_names(excel_names),
           "Output/20220425/GOMeth/BL2ST_Promoter.xlsx",
           format_headers = F)

filepaths[
  enrichment_results_promoter %>% map_lgl(~ filter(.x, FDR < 0.2) %>% nrow %>% `==`(0) %>% `!`)]
enrichment_results_promoter[[
  enrichment_results_promoter %>% map_lgl(~ filter(.x, FDR < 0.2) %>% nrow %>% `==`(0) %>% `!`) %>% which]]

enrichment_results_all <-
  map(filepaths, run_enrichment_test,
      genomicfeaturestest = "ALL")

write_xlsx(enrichment_results_all %>% 
             set_names(excel_names),
           "Output/20220425/GOMeth/BL2ST_All.xlsx",
           format_headers = F)





# long term change ------------------------------------------------------------

filepaths <-
  list.files(path = "Output/20220425/ResultsTables", full.names = T) %>%
  .[grepl("_120_|_180_", .)] %>%
  .[grepl("covarB", .)]
excel_names <-
  sapply(filepaths, convert_filepath_to_outname)

enrichment_results_promoter <-
  map(filepaths, run_enrichment_test,
      genomicfeaturestest = c("TSS1500", "TSS200", "1stExon"))

write_xlsx(enrichment_results_promoter %>% 
             set_names(excel_names),
           "Output/20220425/GOMeth/BL2LT_Promoter.xlsx",
           format_headers = F)

enrichment_results_all <-
  map(filepaths, run_enrichment_test,
      genomicfeaturestest = "ALL")

write_xlsx(enrichment_results_all %>% 
             set_names(excel_names),
           "Output/20220425/GOMeth/BL2LT_All.xlsx",
           format_headers = F)




# baseline ---------------------------------------------------------------------

filepaths <-
  list.files(path = "Output/20220425/ResultsTables", full.names = T) %>%
  .[grepl("BLonly_", .)] %>%
  .[grepl("covarA", .)]
excel_names <-
  sapply(filepaths, convert_filepath_to_outname) %>% gsub("BLonly_", "", .)

enrichment_results_promoter <-
  map(filepaths, run_enrichment_test,
      genomicfeaturestest = c("TSS1500", "TSS200", "1stExon"))

write_xlsx(enrichment_results_promoter %>% 
             set_names(excel_names),
           "Output/20220425/GOMeth/BLonly_Promoter_covarA.xlsx",
           format_headers = F)

enrichment_results_all <-
  map(filepaths, run_enrichment_test,
      genomicfeaturestest = "ALL")

write_xlsx(enrichment_results_all %>% 
             set_names(excel_names),
           "Output/20220425/GOMeth/BLonly_All_covarA.xlsx",
           format_headers = F)






# baseline ---------------------------------------------------------------------

filepaths <-
  list.files(path = "Output/20220425/ResultsTables", full.names = T) %>%
  .[grepl("Trt|Time", .)] %>%
  .[grepl("covarA|covar0", .)]
excel_names <-
  sapply(filepaths, convert_filepath_to_outname)

enrichment_results_promoter <-
  map(filepaths, run_enrichment_test,
      genomicfeaturestest = c("TSS1500", "TSS200", "1stExon"))

write_xlsx(enrichment_results_promoter %>% 
             set_names(excel_names),
           "Output/20220425/GOMeth/TimeTrt_Promoter_covarA.xlsx",
           format_headers = F)

enrichment_results_all <-
  map(filepaths, run_enrichment_test,
      genomicfeaturestest = "ALL")

write_xlsx(enrichment_results_all %>% 
             set_names(excel_names),
           "Output/20220425/GOMeth/TimeTrt_All_covarA.xlsx",
           format_headers = F)






# lmm ---------------------------------------------------------------------

filepaths <-
  list.files(path = "Output/20220425/ResultsTables", full.names = T) %>%
  .[grepl("lmm", .)]
excel_names <-
  sapply(filepaths, convert_filepath_to_outname) %>%
  gsub("lmm_", "", .)

enrichment_results_promoter <-
  map(filepaths, run_enrichment_test,
      genomicfeaturestest = c("TSS1500", "TSS200", "1stExon"))

write_xlsx(enrichment_results_promoter %>% 
             set_names(excel_names),
           "Output/20220425/GOMeth/lmm_Promoter_covarA.xlsx",
           format_headers = F)

enrichment_results_all <-
  map(filepaths, run_enrichment_test,
      genomicfeaturestest = "ALL")

write_xlsx(enrichment_results_all %>% 
             set_names(excel_names),
           "Output/20220425/GOMeth/lmm_All_covarA.xlsx",
           format_headers = F)

