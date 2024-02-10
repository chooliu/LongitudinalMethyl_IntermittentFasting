# ==============================================================================
# D/04b_caret_BLme_to_ShortTerm.R
# model delta 3mo or 6mo outcome ~ baseline methylome
# ==============================================================================


# run interactively ------------------------------------------------------------

screen -S multivar
qlogin -R rusage[mem=12]
cd $HOME/Analyses/DRIFT2_Methylation
module load R/4.1.0_beta
R



# setup ------------------------------------------------------------------------

library(caret)
library(tidyverse)
library(glmnet)

date_export <- "20220425"
dir.create(paste0("Output/", date_export))
dir.create(paste0("Output/", date_export, "/CaretModels/"))



# load & prep data -------------------------------------------------------------

load("Output/20220425/metadata_final.Rdata")
load("Output/20220425/final_methylation_modeling.Rdata")

baseline_zvals <-
  zvals_ruvg_cellcomp[c(T, F), ]

ztransform <-
  function(x) { (x - mean(x, na.rm = T))/sd(x, na.rm = T) }

filepaths <-
  list.files(path = "Output/20220425/ResultsTables") %>%
  .[grepl("Methyl0", .) & grepl("_30_|_60_", .)] %>%
  .[grepl("RUVcellc", .)] %>%
  .[( grepl("covarA", .) & grepl("bmi|weight", .) ) | grepl("covarB", .)]

models_to_run <-
  tibble(target_tsv = filepaths) %>%
  rowwise() %>%
  mutate(out_label =
           str_split(target_tsv, "_") %>%
           unlist %>% .[2:3] %>% paste0(collapse = "_") %>% gsub(".tsv", "", .))




# fxn to loop over features ----------------------------------------------------
# assumes covariate set B (unless bmi or weight, then A)

run_caret_models <- function(target_tsv, out_label) {
  
  # prep data
  
  print(outcomedelta_label)

  outcomedelta_label <-
    gsub("D", "outcomedelta_", out_label)
  print(outcomedelta_label)
  
  resultstable_filtered <-
    read_tsv(paste0("Output/20220425/ResultsTables/", target_tsv)) %>%
    filter(P < 0.001)

  caret_input_dat <-
    baseline_zvals[ , probenames_filtered %in% resultstable_filtered$Probe] %>%
    as_tibble() %>%
    set_names(., probenames_filtered %>% .[`%in%`(., resultstable_filtered$Probe)]) %>% 
    bind_cols(
      metadata_final_subject %>%
        select(!!sym(outcomedelta_label), sex, age, genoPC1, genoPC2, outcome_bmi_0)) %>%
    filter( !is.na(!!sym(outcomedelta_label)) ) %>%
    mutate_if(.predicate = is.numeric, .funs = ztransform)
  
  clinical_only_dat <-
    metadata_final_subject %>%
    select(!!sym(outcomedelta_label),
           Intervention, sex, age, outcome_bmi_0, genoPC1, genoPC2,
           matches("_0$")) %>%
    mutate(Intervention = as.numeric(Intervention == "IMF")) %>%
    filter( !is.na(!!sym(outcomedelta_label)) ) %>%
    mutate_if(.predicate = is.numeric, .funs = ztransform)
  
    filter_subjs_complete3mo <-
      apply(clinical_only_dat, 1, function(x) { is.na(x) %>% any %>% `!` })
    
  print(paste0("N = ", nrow(caret_input_dat)))

  # run models
  
  caret_criteria <-
    trainControl(method = "LOOCV")
  formula_allpred <-
    paste0(outcomedelta_label, " ~ .") %>% as.formula
  
  set.seed(1111)
  print("lasso_fit")
  lasso_fit <-
    train(formula_allpred,
    data = caret_input_dat,
    method = "glmnet",
    trControl = caret_criteria,
    tuneGrid = expand.grid(alpha = 1,
                           lambda = c(seq(1e-6, 0.1, length.out = 50),
                                      seq(0.1, 1, length.out = 50))
  ))
  
  set.seed(1111)
  print("ranger_fit")
  ranger_fit <-
    train(formula_allpred,
          data = caret_input_dat,
          method = "ranger",
          trControl = caret_criteria,
          importance = 'impurity',
          tuneGrid =
            expand.grid(mtry = seq(2, sqrt(ncol(caret_input_dat))*2, length.out = 5) %>% ceiling,
                        splitrule = c("extratrees"), 
                        min.node.size = c(1, 2)))
  
  set.seed(1111)
  print("lasso_fit_filt")
  lasso_fit_filt <-
    train(formula_allpred,
          data = caret_input_dat[filter_subjs_complete3mo, ],
          method = "glmnet",
          trControl = caret_criteria,
          tuneGrid = expand.grid(alpha = 1,
                                 lambda = c(seq(1e-6, 0.1, length.out = 50),
                                            seq(0.1, 1, length.out = 50))
          ))
  
  set.seed(1111)
  print("ranger_fit_filt")
  ranger_fit_filt <-
    train(formula_allpred,
          data = caret_input_dat[filter_subjs_complete3mo, ],
          method = "ranger",
          trControl = caret_criteria,
          importance = 'impurity',
          tuneGrid =
            expand.grid(mtry = seq(2, sqrt(ncol(caret_input_dat))*2, length.out = 5) %>% ceiling,
                        splitrule = c("extratrees"), 
                        min.node.size = c(1, 2)))
  
  set.seed(1111)
  print("lasso_cntl")
  lasso_cntl <-
    train(formula_allpred,
          data = clinical_only_dat[filter_subjs_complete3mo, ],
          method = "glmnet",
          trControl = caret_criteria,
          tuneGrid = expand.grid(alpha = 1,
                                 lambda = c(seq(1e-6, 0.1, length.out = 50),
                                            seq(0.1, 1, length.out = 50))
          ))
  
  set.seed(1111)
  print('ranger_cntl')
  ranger_cntl <-
    train(formula_allpred,
          data = clinical_only_dat[filter_subjs_complete3mo, ],
          method = "ranger",
          trControl = caret_criteria,
          importance = 'impurity',
          tuneGrid =
            expand.grid(mtry = seq(2, sqrt(ncol(clinical_only_dat))*2, length.out = 5) %>% ceiling,
                        splitrule = c("extratrees"), #   Error: Gini splitrule applicable to classification data only.
                        min.node.size = c(1, 2)))
  
  set.seed(1111)
  print('linmodel_cntl')
  linmodel_cntl <-
    train(formula_allpred,
          data = clinical_only_dat[filter_subjs_complete3mo, ],
          method = "lm",
          trControl = caret_criteria)

  # export
  
  save(lasso_fit, ranger_fit, lasso_fit_filt, ranger_fit_filt,
       lasso_cntl, ranger_cntl, linmodel_cntl,
       caret_input_dat,
       file = paste0("./Output/20220425/CaretModels/", out_label, "_RUVcellc_covarAB.Rdata"))
}



# apply helper fxn -------------------------------------------------------------

models_to_run %>%
  arrange(rev(out_label)) %>%
  pmap(run_caret_models)


