# ==============================================================================
# C/05_compile_annot_DMP.R 
# loop through each results file --> get top hits --> annotate --> aggregate
# applies the z-score to beta effect size conversions
# ==============================================================================


# # usually run interactively
# screen -S annotation
# cd $HOME/Analyses/DRIFT2_Methylation
# qlogin -R rusage[mem=16]
# module load R/4.1.0_beta
# R


library(data.table)
library(tidyverse)
library(sesameData)
library(readxl)
library(writexl)
library(GenomicRanges)

dir.create("Output/20220425/CompiledResults/")

filepaths <-
  list.files(
    path = "Output/20220425/ResultsTables",
    pattern = "*tsv",
    full.names = T)
  

load("./Output/20220425/probe_annotations.Rdata")


annotate_results <-
  function(filepath_in, p_threshold = 1e-4) {
    
    results_out <-
      fread(filepath_in) %>%
      filter(bacon_Pz < p_threshold | fdr < 0.10) %>%
      arrange(bacon_Pz) %>%
      transmute(Probe, Effect = Estimate, SE,
                `t-statistic` = t, `p-value` = bacon_Pz, FDR = fdr,
                Signif = if_else(FDR < 0.1, "*", "")) %>%
    right_join(probe_annotations, ., by = "Probe") %>%
    arrange(FDR)
  
  results_out
  }



convert_filepath_to_outname <-
  function(filepath) {
    gsub("Output/20220425/ResultsTables/|RUV_|RUV||\\.tsv", "", filepath) %>%
      gsub("Methyl0", "M0", .) %>%
      gsub("cellcomp", "cellc", .) %>%
      gsub("covar", "set", .)
  }

convert_effect_sizes <-
  function(df, model_type) {
    if (model_type == 1) {
    df_out <-
      df %>%
      mutate(BetaConversion = BetaConversion1,
             EffectBeta = Effect / BetaConversion1) %>%
      select(-BetaConversion1, -BetaConversion2)
    } else {
      df_out <-
        df %>%
        mutate(BetaConversion = BetaConversion2,
               EffectBeta = Effect / BetaConversion2) %>%
        select(-BetaConversion1, -BetaConversion2)
    }
    return(df_out)
  }




# baseline only ----------------------------------------------------------------

targeted_filenames <-
  filepaths %>% .[grepl("BLonly", .)]

results_set_BL <-
  map(targeted_filenames,
      annotate_results)

results_set_BL <-
  map(1:length(results_set_BL),
    function(i) {
      convert_effect_sizes(
        df = results_set_BL[[i]],
        model_type = if_else(grepl("RUVcellc", targeted_filenames[i]), 2, 1)) })

write_xlsx(results_set_BL %>%
             set_names(targeted_filenames %>%
                         convert_filepath_to_outname %>%
                         gsub("BLonly_", "", .)),
           path = "Output/20220425/CompiledResults/DMPs_BLonly_CrossSectional.xlsx",
           format_headers = F
)





# BL Methyl --> Short Term Change in Outcome (3mo, 6mo) ------------------------

targeted_filenames <-
  filepaths %>% .[grepl("Methyl0_", .)] %>% .[grepl("_30_|_60_", .)]

results_set_B0toDelta3mo6mo <-
  map(targeted_filenames,
      annotate_results)

results_set_B0toDelta3mo6mo <-
  map(1:length(results_set_B0toDelta3mo6mo),
         function(i) {
           convert_effect_sizes(
             df = results_set_B0toDelta3mo6mo[[i]],
             model_type = if_else(grepl("RUVcellc", targeted_filenames[i]), 2, 1)) })

write_xlsx(results_set_B0toDelta3mo6mo %>%
             set_names(targeted_filenames %>%
                         convert_filepath_to_outname %>%
                         gsub("Methyl0_", "", .)),
           path = "Output/20220425/CompiledResults/DMPs_BLme_to_ShortTermChange.xlsx",
           format_headers = F
)




# Treatment Effect & Time ------------------------------------------------------

targeted_filenames <-
  filepaths %>% .[grepl("Trt|Time", .)]

results_set_TrtTime <-
  map(targeted_filenames,
      annotate_results)

results_set_TrtTime <-
  map(1:length(results_set_TrtTime),
      function(i) {
        convert_effect_sizes(
          df = results_set_TrtTime[[i]],
          model_type = if_else(grepl("RUVcellc", targeted_filenames[i]), 2, 1)) }) %>%
  map(~ mutate(.x, EffectBeta = Effect * BetaConversion) )

write_xlsx(results_set_TrtTime %>%
             set_names(targeted_filenames %>% convert_filepath_to_outname),
           path = "Output/20220425/CompiledResults/DMPs_Treatment_and_Time.xlsx",
           format_headers = F
)


