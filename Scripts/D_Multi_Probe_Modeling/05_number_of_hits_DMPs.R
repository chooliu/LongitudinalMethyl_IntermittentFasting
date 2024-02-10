# ==============================================================================
# C/05_number_of_hits_DMPs.R 
# loop through each results file --> get top hits --> annotate --> aggregate
# applies the z-score to beta effect size conversions
# ==============================================================================



library(readxl)
target <- "Output/20220425/CompiledResults/DMPs_BLme_to_ShortTermChange.xlsx"
results_shortterm <- map(readxl::excel_sheets(target),
                         ~ read_xlsx(target, sheet = .x) %>% filter(FDR < 0.10))
results_shortterm_p <- map(readxl::excel_sheets(target),
                           ~ read_xlsx(target, sheet = .x) %>% filter(`p-value` < 1e-4))


DMP_numbers <-
  tibble(feature = readxl::excel_sheets(target),
       fdrsignif = results_shortterm %>% map_dbl(nrow),
       marginal = results_shortterm_p %>% map_dbl(nrow)) %>% 
  filter( (grepl("bmi|weight", feature) & grepl("setA", feature) ) |
  (!grepl("bmi|weight", feature) & grepl("setB", feature) ) ) %>%
  # mutate(cellcomp = grepl("cellc", feature)) %>%
  separate(feature, into = c(NA, "outcome", "time", NA, "covar")) %>%
  mutate(outcome = str_sub(outcome, 2, length(outcome))) %>%
  pivot_wider(names_from = "covar", id_cols = c("outcome", "time"), values_from = c("marginal", "fdrsignif")) %>% 
  right_join(tidy_outcome_names, ., by = c("variablename" = "outcome")) %>%
  select(-variablename)




dmr <- "Output/20220425/CompiledResults/DMRs_ShortTermChange.xlsx"

dmrs_p <- map(readxl::excel_sheets(dmr),
              ~ read_xlsx(dmr, sheet = .x) %>%
                filter(`DMR Signif` == "*"))



DMR_numbers <-
  tibble(feature = readxl::excel_sheets(dmr),
       signif = dmrs_p %>% map_dbl(nrow)) %>% 
  filter( (grepl("bmi|weight", feature) & grepl("setA", feature) ) |
            (!grepl("bmi|weight", feature) & grepl("setB", feature) ) ) %>%
  # mutate(cellcomp = grepl("cellc", feature)) %>%
  separate(feature, into = c(NA, "outcome", "time", NA, "covar")) %>%
  mutate(outcome = str_sub(outcome, 2, length(outcome))) %>%
  pivot_wider(names_from = "covar", id_cols = c("outcome", "time"), values_from = c("signif")) %>% 
  right_join(tidy_outcome_names, ., by = c("variablename" = "outcome")) %>%
  select(-variablename)




left_join(DMP_numbers, DMR_numbers) %>% View()
