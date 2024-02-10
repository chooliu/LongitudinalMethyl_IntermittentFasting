# ==============================================================================
# B/07_methyl_qc.R
# check control probes by chip, row, experimental groups of interest
# ==============================================================================




# load data, libraries ---------------------------------------------------------

library(tidyverse)
library(viridis)

load("Output/20220425/metadata_final.Rdata")


# plot features by chip position -----------------------------------------------

plot_by_assay_pos <-
  function(variable, metadat_in = metadata_final_sample) {
  metadat_in %>%
  group_by(EPIC_chip) %>%
  group_split() %>%
  bind_rows() %>%
  mutate(EPIC_row = gsub("C01", "", EPIC_row) %>% as.factor() %>% fct_rev(),
         EPIC_chip = as.factor(EPIC_chip) %>%
           choomisc::set_fctlevels(., c(paste0("P1_", 1:12), paste0("P2_", 1:4)))) %>% 
  ggplot(., aes_string(x = "EPIC_chip", y = "EPIC_row", fill = variable)) +
  geom_tile() + 
  geom_vline(xintercept = c(6.5, 10.5, 16.5)) +
  scale_fill_viridis(option = "D") +
  scale_x_discrete("Plate _ EPIC Chip\n(order in assay plate/order scanned)",
                   expand = c(0, 0)) +
  scale_y_discrete("EPIC Row", expand = c(0, 0)) + 
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90),
        legend.key.width = unit(1.5, "cm"), legend.key.height = unit(0.5, "cm"))
}



# diff intensity patterns by row -----------------------------------------------

# some variation in global methylation levels by chip x row,
# but not consistent patterns (not uniform decrease by row; only some diff chips)
plot_by_assay_pos("mean_beta_cg")
plot_by_assay_pos("mean_beta_ch")


plot_by_assay_pos("mean_intensity")
plot_by_assay_pos("`Column-Scaled Intensity`",
      metadata_final_sample %>%
      group_by(EPIC_chip) %>%
      group_split() %>%
      map(~ mutate(.x, `Column-Scaled Intensity` =
            mean_intensity - mean(mean_intensity))) %>%
            bind_rows())

# some variability in bisulfite conversion metrics
# consider omitting one participant (ID 039), w.r.t. conversion
# metadata_final_sample %>% arrange(-BisulfiteConversion2) %>% head %>% .$IID

plot_by_assay_pos("BisulfiteConversion2")

# other metrics
plot_by_assay_pos("Specificity2")
plot_by_assay_pos("Specificity2Background")

plot_by_assay_pos("NonPolymorphicGreen")
plot_by_assay_pos("NonPolymorphicRed")

plot_by_assay_pos("num_na_cg")
plot_by_assay_pos("num_na_ch")


# plot density plots of sesame qc ----------------------------------------------

plot_sesame_qc_metrics <-
  function(var_to_plot) {
    plot <-
      ggplot(metadata_final_subject,
             aes_string(x = var_to_plot)) +
      geom_density() +
      theme_ridges() +
      scale_color_manual(values = c("black", "red"))
    plot
  }


plot_sesame_qc_metrics("num_na_cg")
plot_sesame_qc_metrics("num_na_ch")

plot_sesame_qc_metrics("mean_intensity")

plot_sesame_qc_metrics("frac_meth_cg")

plot_sesame_qc_metrics("Specificity2")
plot_sesame_qc_metrics("Specificity2Background")



# inferred participant features via sesame -------------------------------------

# age highly consistent
ggplot(metadata_final_sample,
       aes(InferredAge, age)) +
  geom_point() +
  theme_few()

cor(metadata_final_sample$age,
    metadata_final_sample$InferredAge)

# sex highly consistent
table(metadata_final_sample$InferredSex,
      metadata_final_sample$sex)

table(metadata_final_sample$InferredKarotype,
      metadata_final_sample$sex)

# ancestry
table(metadata_final_sample$ethnicity,
      metadata_final_sample$InferredEthnicity)





# double check features by chip ------------------------------------------------
# seem well-randomized

table1_by_batch <- 
  function(groupvar) {
    metadata_final_subject %>%
      group_by_at(groupvar) %>%
      summarize(n = n(),
                `Self-Reported Race (white)` = summaryCountPercent(race, 5),
                `GenoPC1` = summaryMeanSD(genoPC1),
                `GenoPC2` = summaryMeanSD(genoPC2),
                `Age (years)` = summaryMeanSD(age),
                `Sex (female)` = summaryCountPercent(sex, "0"),
                `Baseline Weight (lb)` = summaryMeanSD(outcome_weight_0),
                `Delta Weight (lb)` = summaryMeanSD(outcomedelta_weight_30),
                `Intervention (IMF)` = summaryCountPercent(Intervention, "IMF"))
}

# EPIC_WellPairs (row) identified from the PCA as possible tech effect
# no obvious difference in demographics by wellpair
table1_by_batch("EPIC_WellPairs")

# also seem balanced by sex, age, intervention, baseline weight, delta weight
table1_by_batch("EPIC_Date_Scanned")
table1_by_batch("EPIC_chip") 


