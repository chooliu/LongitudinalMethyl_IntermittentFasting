# ==============================================================================
# D/08_methylglm_gsea.R
# pathway enrichment via methylglm
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
library(methylGSA)

dir.create("Output/20220425/methylGSA/")
load("Output/20220425/probe_annotations.Rdata")

probe_annot_methylGSAfmt <-
  probe_annotations %>%
  select(Probe, Gene) %>%
  separate_rows(Gene, sep = ";") %>%
  arrange(grepl("^MIR|^LINC|-", Gene, ignore.case = T)) %>% # deprioritize ncRNAs, all else even
  filter(!duplicated(Probe)) %>%
  mutate_all(as.character) %>%
  as.matrix() %>%
  methylGSA::prepareAnnot()



run_methylglm <- function(filepath_in) {
  
  results_table <- fread(filepath_in)
  
  GSEA_Results <-
    methylglm(results_table$bacon_Pz %>% set_names(results_table$Probe) %>% .[!is.na(.)],
              array.type = "EPIC", GS.type = "Reactome",
              FullAnnot = probe_annot_methylGSAfmt,
              minsize = 25, maxsize = 500
    )
  
  GSEA_Results
}


convert_filepath_to_outname <-
  function(filepath) {
    basename(filepath) %>%
      gsub(".tsv", "", .) %>%
      gsub("Methyl0", "M0", .) %>%
      gsub("cellcomp", "cellc", .) %>%
      gsub("covar", "set", .) %>%
      gsub("RUV", "", .) 
  }


# short term

filepaths <-
  list.files(path = "Output/20220425/ResultsTables", full.names = T) %>%
  .[grepl("_30_|_60_", .)] %>%
  .[grepl("covarB", .) | grepl("covarA", .) & grepl("bmi|weight", .)] %>%
  .[grepl("RUVcellc", .)]

results_methylglm <-
  map(filepaths, run_methylglm)

write_xlsx(results_methylglm %>%
             set_names(filepaths %>% convert_filepath_to_outname),
           path = "Output/20220425/methylGSA/Reactome_BLtoST.xlsx",
           format_headers = F
)



# long term --------------------------------------------------------------------

filepaths <-
  list.files(path = "Output/20220425/ResultsTables", full.names = T) %>%
  .[grepl("_120_|_180_", .)] %>%
  .[grepl("covarB", .) | grepl("covarA", .) & grepl("bmi|weight", .)] %>%
  .[grepl("RUVcellc", .)]

results_methylglm <-
  map(filepaths, run_methylglm)

write_xlsx(results_methylglm %>%
             set_names(filepaths %>% convert_filepath_to_outname),
           path = "Output/20220425/methylGSA/Reactome_BLtoLT.xlsx",
           format_headers = F
)



# long term --------------------------------------------------------------------

filepaths <-
  list.files(path = "Output/20220425/ResultsTables", full.names = T) %>%
  .[grepl("_120_|_180_", .)] %>%
  .[grepl("covarB", .) | grepl("covarA", .) & grepl("bmi|weight", .)] %>%
  .[grepl("RUVcellc", .)]

results_methylglm <-
  map(filepaths, run_methylglm)

write_xlsx(results_methylglm %>%
             set_names(filepaths %>% convert_filepath_to_outname),
           path = "Output/20220425/methylGSA/Reactome_BLtoLT.xlsx",
           format_headers = F
)




# trt/time ---------------------------------------------------------------------

filepaths <-
  list.files(path = "Output/20220425/ResultsTables", full.names = T) %>%
  .[grepl("Trt|Time", .)] %>%
  .[grepl("RUVcellc", .)]

results_methylglm <-
  map(filepaths, run_methylglm)

write_xlsx(results_methylglm %>%
             set_names(filepaths %>% convert_filepath_to_outname),
           path = "Output/20220425/methylGSA/Reactome_TrtTime.xlsx",
           format_headers = F
)




