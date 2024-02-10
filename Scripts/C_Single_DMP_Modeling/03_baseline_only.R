# ==============================================================================
# C/03_baseline_only.R 
# examine cross-sectional baseline methylation <--> cardiometabolic features
# ==============================================================================




# libraries & input data -------------------------------------------------------

library(tidyverse)
library(bacon)

date_export <- "20220425"
dir.create(paste0("Output/", date_export))
dir.create(paste0("Output/", date_export, "/ResultsTables/"))
dir.create(paste0("Output/", date_export, "/BaconParams/"))

load("Output/20220425/metadata_final.Rdata")
load("Output/20220425/final_methylation_modeling.Rdata")

# modelindex from R CMD BATCH/bsub (see 00_index.sh file)
# (to run a subset of models, in one of ten batches)
modelindex <- commandArgs(trailingOnly = T) %>% as.numeric
nbatch <- 10




# core modeling fxn ------------------------------------------------------------

run_linear_model <-
  function(formula, # linear model formula
           metadata_table = metadata_final_subject, # sample or subject-wise
           outputname, # results table export name
           df_in, # degrees of freedom for t --> z-statistic
           index_for_effect = 2, # get 2nd regression coef by default
           seed_in = 12345, # set seed for bacon
           save_in_global_env = F, # optional export to R and
           max_probes = 1e6) { # change to small probe # for troubleshooting
    
    print(paste0("running... '", outputname, "'"))
    
    # running linear models
    max_probes <-
      min(max_probes, ncol(zvals_ruvg))
    
    # extract effects
    results <-
      sapply(1:max_probes,
             function(y) {
               lm(
                 formula %>% as.formula, 
                 metadata_table) %>%
                 summary %>% coef %>% .[index_for_effect, ]
             })
    
    # summarize each probe's effect size
    results <-
      results %>%
      t %>% as_tibble() %>%
      set_names(c("Estimate", "SE", "t", "P")) %>%
      bind_cols(Probe = probenames_filtered[1:max_probes], .) %>%
      arrange(P) %>%
      rowwise() %>%
      mutate(z_stat = qnorm( pt( t, df = df_in) )) # *
    
    # running bacon
    # [*] rare cases with very small p-values;
    # for t > 11.6, df ~ 55,  returns z-score = Inf
    # t <= 11.6 returns z-score = 8.21, so replace with this value
    # (likely an underestimate / erring on being conservative)
    results$z_stat[!is.finite(results$z_stat)] <- 8.21
    filter_P_ok <-
      !is.na(results$P) & is.finite(results$z_stat)
    
    # run bacon for epigenome-wide bias and inflation correction
    set.seed(seed_in)
    bacon_obj <- bacon(results$z_stat[filter_P_ok])
    results$bacon_Pz <- NA_real_
    results$bacon_Pz[filter_P_ok] <- bacon_obj %>% pval %>% as.numeric
    results$fdr <- p.adjust(results$bacon_Pz, "fdr")
    
    # exporting bacon parameters & results table
    # results table: Probe, Estimate, SE, t, P, z_stat, bacon_Pz, fdr
    write_tsv(
      c(estimates(bacon_obj)) %>% enframe %>% set_names(c("Param", "Value")),
      paste0("./Output/", date_export, "/BaconParams/", outputname, ".txt"))
    
    write_tsv(
      results,
      file = paste0("./Output/", date_export, "/ResultsTables/", outputname, ".tsv"))
    
    # optionally export results into R global environment
    # for viz / pipeline troubleshooting
    if (save_in_global_env) {
      assign(x = paste0("results_", outputname),
             value = results, envir = .GlobalEnv)
    }
    
    print(paste0("finished running: 'results_", outputname, "'"))
    
  }




# generate model LHS/RHS combos ------------------------------------------------

# RHS covariates
covariate_set_names <-
  tibble(
    covariate_set = paste0("covar", c("A", "B", "D")),
    RHS = c(" + sex + age + genoPC1 + genoPC2",
            " + sex + age + genoPC1 + genoPC2 + outcome_bmi_0",
            " + sex + age + genoPC1 + genoPC2 + outcome_bmi_0 + CellTypePC1 + CellTypePC2"),
    df = 62 - 1 - c(4, 5, 7)) # subtract df for intercept & covariates

# LHS labels
# mvals_ruvg[c(T, F), ] = baseline values, to baseline cardiometabolic
outcome_set_names <-
  metadata_final_subject %>% names %>% .[grepl("outcome", .) & grepl("_0", .)] %>%
  sapply(function(x) {
    c(prefix = paste0("BLonly_", gsub("outcome_|_0", "", x), "_RUV"),
      LHS = paste0(x, " ~ zvals_ruvg[c(T, F), y, drop = T]"))
  }) %>%
  t %>%
  as_tibble() 




# compile LHS & RHS ------------------------------------------------------------
# below generates each combination of LHS & RHS

conditions_to_run <-
  expand_grid(LHS = outcome_set_names$LHS,
              RHS = covariate_set_names$RHS) %>%
  left_join(covariate_set_names) %>%
  left_join(outcome_set_names) %>%
  bind_rows(., mutate(., prefix = paste0(prefix, "cellc"), # run +/- cell comp
                      LHS = gsub("zvals_ruvg", "zvals_ruvg_cellcomp", LHS))) %>%
  filter( !(covariate_set != "covarA" & grepl("bmi|weight", prefix)) ) %>%  # [1] baseline BMI in model 2x
  filter(! (grepl("cellcomp", LHS) & covariate_set == "covarD") ) %>% # [2] only correct for cell comp in the z-scores w/o existing cc adjustment
  transmute(formula = paste0(LHS, RHS),
            outputname = paste0(prefix, "_", covariate_set),
            df_in = df)




# run a subset of models (1 of 10 batches) --------------------------------------
# 88 conditions, break into sections to partially parallelize

filter_models_to_run <-
  rep(c(1:nbatch), each = ceiling(nrow(conditions_to_run)/nbatch)) %>%
  .[1:nrow(conditions_to_run)] %>%
  `==`(modelindex)

pmap(conditions_to_run[filter_models_to_run, ],
     run_linear_model,
     metadata_table = metadata_final_sample %>% filter(Timepoint == "BL"))



