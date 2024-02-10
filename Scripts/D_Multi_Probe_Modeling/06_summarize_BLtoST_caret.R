# ==============================================================================
# D/06_summarize_BLtoST_caret.R
# loop through the BL --> short-term outcome models, summarize performance
# ==============================================================================



# fast -- run interactively ----------------------------------------------------

screen -S multivar
qlogin -R rusage[mem=12]
cd $HOME/Analyses/DRIFT2_Methylation
module load R/4.1.0_beta
R


# libraries --------------------------------------------------------------------

library(caret)
library(tidyverse)
library(glmnet)
library(ranger)
library(pROC)
library(writexl)

load("./Output/20220425/probe_annotations.Rdata")

filepaths_modelout <-
  list.files("Output/20220425/CaretModels/", pattern = "*.Rdata",
             full.names = T) %>%
  .[grepl("_30_|_60_", .)]  %>% 
  .[grepl("covarAB", .)]



# helper fxn to extract model info ---------------------------------------------

format_caret_results <-
  function(filepath_in) {
    
    load(filepath_in)
    print(filepath_in)

    # summarizing variable importance
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
      left_join(probe_annotations, by = "Probe") %>%
      arrange(LassoRank, RangerRank)
    
    
    # final prediction metrics
    if (lasso_fit$modelType == "Regression") {
      
      across_model_performance <-
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
        LM = linmodel_cntl$results %>% as.tibble) %>%
        map_df( ~ select(.x, RMSE, Rsquared, MAE)) %>%
        bind_cols(Method = c("LASSO", "RF",  "LASSO_filt", "RF_filt", "LASSO_cntl", "RF_cntl", "LM_cntl"), .)
      
    } else {
      return(NA)
    }
    
    # return variable importance & performance
    list(VarImp = final_caret_varimp,
         Perform = across_model_performance)
    
}



# apply helper fxn -------------------------------------------------------------
# export variable importance & performance (r2, RSME, MAE)

all_caret_results <-
  filepaths_modelout %>% map(format_caret_results)

excel_names <-
  filepaths_modelout %>% basename(.) %>%
  gsub(".Rdata", "", .) %>%
  gsub("RUVcellc_", "Rcc_", .) 

all_caret_results %>% map(~ .[["VarImp"]]) %>% set_names(excel_names) %>%
  map(~ mutate(.x, BetaConversion = BetaConversion2) %>%
        select(-BetaConversion1, -BetaConversion2)) %>%
  write_xlsx(., path = "./Output/2023/BLtoST_Multivariate_FeatSelection_covarAB.xlsx", format_headers = F)

map(1:length(excel_names),
    function(i) {
      all_caret_results[[i]][["Perform"]] <<-
        all_caret_results[[i]][["Perform"]] %>%
        bind_cols(Feature = excel_names[i], .)
    }) 


all_caret_results %>% map(~ .[["Perform"]]) %>%
  map(~ pivot_longer(.x, cols = 3:5)) %>%
  map_dfr(~ pivot_wider(.x, names_from = c("Method", "name"), values_from = "value")) %>%
  write_xlsx(., path = "./Output/2023/BLtoST_Performance_covarAB.xlsx", format_headers = F)



# troubleshooting note about cntl r2 for certain models: -----------------------
# some can yield impossible r2 = 1 because LASSO selects an intercept-only model

# lasso_cntl_vars <-
#   coef(lasso_cntl$finalModel, lasso_cntl$bestTune$lambda) %>%
#   as.matrix # sparse matrix --> mat --> tibble
