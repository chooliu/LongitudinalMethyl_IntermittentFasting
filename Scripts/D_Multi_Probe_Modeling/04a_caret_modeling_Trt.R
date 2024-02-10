# ==============================================================================
# D/04a_caret_modeling_Trt.R
# model Trt group ~ delta 3mo methylation
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
library(pROC)

date_export <- "20220425"
dir.create(paste0("Output/", date_export))
dir.create(paste0("Output/", date_export, "/CaretModels/"))

load("Output/20220425/metadata_final.Rdata")
load("Output/20220425/probe_annotations.Rdata")
load("Output/20220425/final_methylation_modeling.Rdata")




# load & prep data -------------------------------------------------------------

ztransform <-
  function(x) { (x - mean(x, na.rm = T))/sd(x, na.rm = T) }

delta_zvals <-
  zvals_ruvg_cellcomp[c(F, T), ] -
  zvals_ruvg_cellcomp[c(T, F), ]

resultstable_filtered <-
  read_tsv("./Output/20220425/ResultsTables/Trt_DMethyl_RUVcellc_covarB.tsv") %>%
  filter(P < 0.001)

caret_input_dat <-
  delta_zvals[ , probenames_filtered %in% resultstable_filtered$Probe] %>%
  as_tibble() %>%
  set_names(., probenames_filtered %>% .[`%in%`(., resultstable_filtered$Probe)]) %>% 
  bind_cols(metadata_final_subject %>%
              select(Intervention, sex, age, genoPC1, genoPC2, outcome_bmi_0),
            .) %>%
  mutate(Intervention = Intervention %>% as.factor) %>%
  mutate_if(.predicate = is.numeric, .funs = ztransform)

clinical_only_dat <-
  metadata_final_subject %>%
  select(Intervention, sex, age, outcome_bmi_0, genoPC1, genoPC2,
         matches("_30$")) %>%
  mutate(Intervention = Intervention %>% as.factor) %>%
  mutate_if(.predicate = is.numeric, .funs = ztransform)


filter_subjs_complete3mo <-
  apply(clinical_only_dat, 1, function(x) { is.na(x) %>% any %>% `!` })




# run models --> export --------------------------------------------------------
  
caret_criteria <-
  trainControl(method = "LOOCV",
               summaryFunction = twoClassSummary,
               classProbs = TRUE
  )
formula_allpred <- "Intervention ~ ." %>% as.formula

set.seed(1111)
lasso_fit <-
  train(formula_allpred,
  data = caret_input_dat,
  method = "glmnet",
  verbose = FALSE,
  metric = "ROC",
  trControl = caret_criteria,
  family = "binomial",
  tuneGrid = expand.grid(
    alpha = 1, 
    lambda = c(seq(1e-6, 0.1, length.out = 50), seq(0.1, 1, length.out = 50))
  ))

set.seed(1111)
ranger_fit <-
  train(formula_allpred,
        data = caret_input_dat,
        method = "ranger",
        verbose = FALSE,
        metric = "ROC",
        num.trees = 200,
        trControl = caret_criteria,
        importance = 'impurity',
        tuneGrid = 
        expand.grid(mtry = seq(2, sqrt(ncol(caret_input_dat))*2, length.out = 5) %>% ceiling,
                    splitrule = c("extratrees"),
                    min.node.size = c(1, 2)))

set.seed(1111)
lasso_fit_filt <-
  train(formula_allpred,
        data = caret_input_dat[filter_subjs_complete3mo, ],
        method = "glmnet",
        verbose = FALSE,
        metric = "ROC",
        trControl = caret_criteria,
        family = "binomial",
        tuneGrid = expand.grid(
          alpha = 1, 
          lambda = c(seq(1e-6, 0.1, length.out = 50), seq(0.1, 1, length.out = 50))
        ))

set.seed(1111)
ranger_fit_filt <-
  train(formula_allpred,
        data = caret_input_dat[filter_subjs_complete3mo, ],
        method = "ranger",
        verbose = FALSE,
        metric = "ROC",
        num.trees = 200,
        trControl = caret_criteria,
        importance = 'impurity',
        tuneGrid = 
          expand.grid(mtry = seq(2, sqrt(ncol(caret_input_dat))*2, length.out = 5) %>% ceiling,
                      splitrule = c("extratrees"),
                      min.node.size = c(1, 2)))


set.seed(1111)
lasso_cntl <-
  train(formula_allpred,
        data = clinical_only_dat[filter_subjs_complete3mo, ],
        method = "glmnet",
        trControl = caret_criteria,
        family = "binomial",
        tuneGrid = expand.grid(
          alpha = 1, 
          lambda = c(seq(1e-6, 0.1, length.out = 50), seq(0.1, 1, length.out = 50))
        ))

set.seed(1111)
ranger_cntl <-
  train(formula_allpred,
        data = clinical_only_dat[filter_subjs_complete3mo, ],
        method = "ranger",
        trControl = caret_criteria,
        importance = 'impurity',
        tuneGrid = 
          expand.grid(mtry = seq(2, sqrt(ncol(clinical_only_dat))*2, length.out = 5) %>% ceiling,
                      splitrule = c("extratrees"), # * 
                      min.node.size = c(1, 2)))

set.seed(1111)
cntl_logistic_regression <-
  train(formula_allpred,
        data = clinical_only_dat[filter_subjs_complete3mo, ],
        method = "glm",
        family = "binomial",
        trControl = caret_criteria)

save(lasso_fit, ranger_fit, lasso_fit_filt, ranger_fit_filt,
     lasso_cntl, ranger_cntl, cntl_logistic_regression,
     caret_input_dat,
     file = "./Output/20220425/CaretModels/Trt_DMethyl_RUVcellc_covarB.Rdata")



# summarizing variable importance ----------------------------------------------

load("./Output/20220425/CaretModels/Trt_DMethyl_RUVcellc_covarB.Rdata")
library(writexl)

lasso_df <-
  coef(lasso_fit$finalModel, lasso_fit$bestTune$lambda) %>%
  as.matrix %>% # sparse matrix --> mat --> tibble
  as_tibble(rownames = "Probe") %>% 
  filter(s1 != 0) %>%
  filter(Probe != "(Intercept)") %>%
  arrange(-abs(s1)) %>% # arrange(Probe != "(Intercept)", -abs(`1`)) %>%
  bind_cols(., LassoRank = 1:nrow(.))
names(lasso_df)[2] <- "LassoCoef"

ranger_df <-
  ranger_fit %>% varImp(scale = F) %>% .$importance %>%
  filter(Overall >= quantile(.$Overall, prob = 0.75)) %>%
  as_tibble(rownames = "Probe") %>%
  arrange(-Overall) %>%
  bind_cols(., RangerRank = 1:nrow(.))
names(ranger_df)[2] <- "Ranger Gini"

final_caret_varimp <-
  full_join(lasso_df, ranger_df, by = "Probe") %>%
  left_join(probe_annotations %>% mutate(BetaConversion = BetaConversion2) %>%
              select(-BetaConversion1, -BetaConversion2), by = "Probe") %>%
  arrange(LassoRank, RangerRank)

write_xlsx(final_caret_varimp,
           path = "./Output/20220425/Trt_Multivariate_FeatSelection_covarB.xlsx",
           format_headers = F)



# summarizing model performance ------------------------------------------------

list(
  LASSO = lasso_fit$results %>% as_tibble() %>% filter(lambda == lasso_fit$best$lambda),
  Ranger = ranger_fit$results %>% as_tibble() %>%
    filter(mtry == ranger_fit$best$mtry & min.node.size == ranger_fit$best$min.node.size),
  LASSO_filt = lasso_fit_filt$results %>% as_tibble() %>% filter(lambda == lasso_fit$best$lambda),
  Ranger_filt = ranger_fit_filt$results %>% as_tibble() %>%
    filter(mtry == ranger_fit$best$mtry & min.node.size == ranger_fit$best$min.node.size),
  Lasso_Cntl = lasso_cntl$results %>% as_tibble() %>% filter(lambda == lasso_cntl$best$lambda),
  Ranger_Cntl = ranger_cntl$results %>% as_tibble() %>%
    filter(mtry == ranger_cntl$best$mtry & min.node.size == ranger_cntl$best$min.node.size),
  Logistic_Cntl = cntl_logistic_regression$results %>% as.tibble) %>%
  map_df( ~ select(.x, ROC, Sens, Spec)) %>%
  bind_cols(Model = c("LASSO", "RF", "LASSOfilt", "RFfilt", "LASSO_cntl", "RF_cntl", "Logistic_cntl"), .) %>%
  write_xlsx(., path = "./Output/20220425/Trt_Performance_covarB.xlsx", format_headers = F)
