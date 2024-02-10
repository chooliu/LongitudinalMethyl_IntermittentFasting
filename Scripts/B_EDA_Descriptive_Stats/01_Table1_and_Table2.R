# ==============================================================================
# B/01_Table1_and_Table2
# check demographics and outcomes
# ==============================================================================




# load data, libraries ---------------------------------------------------------
# note: these are needed for all B/ scripts

library(tidyverse)
library(ggthemes)
library(ggridges)
library(cowplot)
library(RColorBrewer)

load("./Output/20220425/metadata_final.Rdata")



# create tidy variable names <--> clean name conversion ------------------------

tidy_outcome_names <-
  metadata_final_subject %>%
  select(contains("outcomedelta") & ends_with("0")) %>%
  names %>%
  map_chr(~ str_split(.x, "_", simplify = T) %>% .[[2]]) %>%
  sort %>% unique %>%
  tibble(variablename = .)
tidy_outcome_names$Feature <-
  c("BMI (kg/m^2)", "Total Cholesterol (mg/dl)", "Cortisol (ug/dL)",
    "C-Reactive Protein (mg/L)", "Diastolic BP (mm Hg)",
    "Trunk Fat (%)", "Whole Body Fat (%)", "Total Trunk Fat (lb)",
    "Whole Body Fat (lb)", "Free Fatty Acid (uEq/L)",
    "Ghrelin (pg/mL)", "Glucose (mg/dL)", "High Density Lipoprotein (mg/dL)",
    "Hemoglobin A1c (%)", "Insulin (uIU/mL)",
    "Low-Density Lipoprotein (mg/dl)", "Leptin (ng/nL)",
    "Pancreatic Peptide YY (pg/mL)", "Systolic BP (mm Hg)",
    "Triglycerides (mg/dL)", "Waist Circumference (cm)", "Weight (lb)")




# Table 1 & summary stat functions ---------------------------------------------


summaryCountPercent <-
  function(x, values, count_NAs = T,
         digits = 1,
         NA_percent_in_brackets = F,
         NA_count_in_brackets = F,
         NA_count_in_brackets_alt = T,
         fuzzy = F, inverse = F) {
  
  n_NA <- sum(is.na(x))
  ratio_miss <- n_NA/length(x)*100
  
  if (n_NA != 0) {
    warning(paste0(n_NA,
                   " values are missing/NA in input vector"))
  }
  
  if (n_NA == length(x)) {
    warning("all values missing!")
    return("-")
  }
  
  # exact match (default)
  n <- sum(x %in% values, na.rm = T)
  tot <- ifelse(count_NAs, length(x), length(x) - n_NA)
  
  # fuzzy matching
  if (fuzzy) {
    search_term <-
      paste(values, collapse = "|")
    n <- grepl(search_term, x, ignore.case = T) %>% sum
  }
  
  # inverse
  if (inverse) {
    n <- tot - n
  }
  
  output <- paste0(n, " (", formatC(n / tot * 100, format = "f", digits = digits),
                   "%)")# (", n, "/", tot, ")")
  
  if (NA_count_in_brackets & n_NA != 0) {
    output <-
      paste0(output, " [", n_NA, "]")
  }
  
  if (NA_percent_in_brackets & n_NA != 0) {
    output <-
      ratio_miss %>%
      formatC(., format = "f", digits = digits) %>%
      paste0(output, " [", ., "%]")
  }
  
  return(output)
}

summaryMeanSD <-
  function(x, digits = 1, na.rm = F,
         NA_percent_in_brackets = F,
         NA_count_in_brackets = F,
         NA_count_in_brackets_alt = T) {
  
  if ( !(typeof(x) %in% c("double", "integer")) ) {
    warning("input vector not numeric")
    x <- as.numeric(x)
  }
  tot <- length(x)
  n_NA <- sum(is.na(x))
  ratio_miss <- n_NA/length(x)*100
  
  if (n_NA != 0) {
    warning(paste0(n_NA,
                   " values are missing/NA in input vector"))
  }
  
  if (n_NA == length(x)) {
    warning("all values missing!")
    return("-")
  }
  
  if (na.rm == T) {
    x <- x[!is.na(x)]
  }
  
  output <- c(mean(x), sd(x))
  output <- formatC(output, format = "f", digits = digits)
  
  output <- paste0( output[1], " (", output[2], ")")
  
  if (NA_count_in_brackets & n_NA != 0) {
    output <-
      paste0(output, " [", n_NA, "]")
  }
  
  if (NA_count_in_brackets_alt & !NA_count_in_brackets & n_NA != 0) {
    output <-
      paste0(output, " [n = ", tot - n_NA, "]")
  }
  
  if (NA_percent_in_brackets & n_NA != 0) {
    output <-
      ratio_miss %>%
      formatC(., format = "f", digits = digits) %>%
      paste0(output, " [", ., "%]")
  }
  
  return(output)
  
}


(
Table1_A <-
  metadata_final_subject %>%
    bind_rows(metadata_final_subject %>% mutate(Intervention = "All")) %>%
  group_by(Intervention) %>%
  summarize(n = n(),
            `Self-Reported Race (white)` = summaryCountPercent(race, 5),
            `Self-Reported Ethnicity (non-hispanic)` = summaryCountPercent(ethnicity, 2),
            `Age (years)` = summaryMeanSD(age),
            `Sex (female)` = summaryCountPercent(sex, "0"),
            `Baseline BMI (kg/m2)` = summaryMeanSD(outcome_bmi_0)
            ) %>%
  mutate_all(.funs = as.character) %>%
  pivot_longer(cols = 2:ncol(.)) %>%
  pivot_wider(id_cols = name,
              names_from = Intervention)
) %>%
  transmute(` ` = name, DCR, IMF, All)


run_lin_model_on_change <-
  # function(y) { lm(y ~ 1 + sex + age,
  #                  metadata_final_subject) %>%
  #     summary %>% coef %>% .[1, 4] %>% format(format = "fg", digits = 2) }
    function(y) { t.test(y)$p.value %>%
        format(format = "fg", digits = 2) }



# create tidy variable names <--> clean name conversion ------------------------

(
  Table1_CHANGES <-
    metadata_final_subject %>%
    # select(contains("outcomedelta") & ends_with("0")) %>% 
    select((contains("outcome") & ends_with("0")) | contains("outcomedelta")) %>%
    summarize_all(list(N = function(x) { is.na(x) %>% `!` %>% sum %>% as.character() },
                       MeanSD = function(x) summaryMeanSD(x, na.rm = T, NA_count_in_brackets_alt = F),
                       CV = function(x) { formatC(sd(x, na.rm = T) / mean(x, na.rm = T), digits = 2) },
                       Pval = run_lin_model_on_change)) %>%
    # mutate_all(.funs = as.character) %>%
    pivot_longer(cols = 1:ncol(.)) %>%
    mutate(name = gsub("outcomedelta_|outcome_", "", name)) %>%
    separate(col = name, into = c("Feature", "TimePoint", "Test"), sep = "_") %>%
    pivot_wider(id_cols = c(Feature, TimePoint), names_from = Test) %>%
    mutate(TimePoint = gsub("0$", "mo", TimePoint) %>% paste0("Δ", .) %>%
             if_else(`==`(., "Δmo"), "BL", .))
) %>%
  filter(TimePoint %in% c("BL", "Δ3mo", "Δ6mo")) %>%
  right_join(tidy_outcome_names, ., by = c("variablename" = "Feature")) %>%
  select(-variablename) %>%
  pivot_wider(id_cols = Feature, names_from = TimePoint,
              values_from = c("N", "MeanSD", "CV", "Pval")) %>%
  select(Feature, contains("BL"), contains("Δ3mo"), contains("Δ6mo")) %>% 
  arrange(Feature) %>%
  View()



# check short-term outcome differences by intervention -------------------------
# (intervention --> outcome was blinded during the study/analysis,
# checking post-hoc among our subset of participants)

run_lin_model_on_change <-
  function(y) { lm(y ~ 1 + Intervention, # + sex + age + outcome_bmi_0,
                   metadata_final_subject) %>%
      summary %>% coef %>% .[2, 4] } #%>% format(format = "fg", digits = 2) }


# intervention effect, three p-values < 0.05
intervention_pvals <- 
  metadata_final_subject %>%
  select(contains("outcomedelta") & ( ends_with("30") | ends_with("60") )) %>%
  summarize_all(run_lin_model_on_change) %>% 
  pivot_longer(cols = 1:ncol(.)) %>%
  set_names(c("outcomedelta", "pval")) %>%
  arrange(pval)

# none significant after FDR correction
intervention_pvals$pval %>% p.adjust("fdr") %>% head

# check effect
ggplot(metadata_final_subject,
       aes_string("Intervention",
                  intervention_pvals$outcomedelta[1])) +
  geom_boxplot() + 
  theme_few()




# plot distribution of change values -------------------------------------------
# waterfall plot & change-in-outcomes

make_waterfall <-
  function(variable_of_interest,
           yaxislabel = variable_of_interest,
           metadatain = metadata_final_subject) {

  tmpplot <-
    metadatain %>%
    ungroup %>%
    mutate( yaxis := !!sym(variable_of_interest)) %>%
    arrange(-yaxis) %>%
    mutate(roworder = 1:nrow(.),
              outcome_is_pos = yaxis > 0)
      
  maxval <- tmpplot$yaxis %>% abs %>% max(na.rm = T) * 1.1
    ggplot(data = tmpplot,
           aes(x = roworder,
               y = yaxis,
               fill = outcome_is_pos)) +
    geom_hline(yintercept = 0, alpha = 0.8, lty = 1) +
    geom_bar(stat = "identity", width = 0.75, alpha = 0.8) +
    theme_few() +
    scale_x_continuous(name = "", labels = NULL, breaks = NULL) +
      scale_fill_manual(values = c(`TRUE` = "#d73027", `FALSE` = "#4575b4")) +
    scale_y_continuous(yaxislabel, limits = c(-maxval, maxval),
                       breaks = scales::pretty_breaks(n = 8)) +
    theme(legend.position = "none")
  }


make_waterfall("outcomedelta_weight_30", "Δ Weight (lb, 3mo)")
make_waterfall("outcomedelta_tg_60", "Δ Triglycerides (mg/dl, 6mo)")




# plot distribution of changes as density plot ---------------------------------

plot_grid(
  ggplot(metadata_final_subject, aes(x = outcomedelta_tg_60, y = 1)) +
    geom_density_ridges(
      jittered_points = TRUE,
      position = position_points_jitter(width = 0.05, height = 0),
      point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7,
    ) +
    theme_ridges() +
    xlab("Δ Triglycerides") + ylab("Density"),
  
  
  ggplot(metadata_final_subject, aes(x = outcome_tg_6, y = 1)) +
    geom_density_ridges(
      jittered_points = TRUE,
      position = position_points_jitter(width = 0.05, height = 0),
      point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7,
    ) +
    theme_ridges() +
    xlab("6-Month Triglycerides") + ylab("Density"),
  
  
  ggplot(metadata_final_subject, aes(x = outcomedelta_insulin_60, y = 1)) +
    geom_density_ridges(
      jittered_points = TRUE,
      position = position_points_jitter(width = 0.05, height = 0),
      point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7,
    ) +
    theme_ridges() +
    xlab("Δ Insulin") + ylab("Density")
  
  ,
  ggplot(metadata_final_subject, aes(x = outcome_insulin_6, y = 1)) +
    geom_density_ridges(
      jittered_points = TRUE,
      position = position_points_jitter(width = 0.05, height = 0),
      point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7,
    ) +
    theme_ridges() +
      xlab("6-Month Insulin") + ylab("Density"))




