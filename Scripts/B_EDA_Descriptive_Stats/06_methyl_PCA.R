# ==============================================================================
# B/06_methyl_PCA.R
# global patterns in DNAme
# ==============================================================================




# load data, libraries ---------------------------------------------------------

library(tidyverse)
library(ggthemes)
library(cowplot)
library(viridis)

load("./Output/20220421/probe_annotations.Rdata")
load("./Output/20220421/final_methylation_modeling.Rdata")
load("./Output/20220421/metadata_final.Rdata")

# load("Output/20220425/metadata_final.Rdata")
# load("ModelInput/final_methylation_full.Rdata")

pca_probes_include <-
  probe_annotations %>%
  filter(Chr %in% c(1:22)) %>%
  transmute(Probe,
            IQR = Beta_Q75 - Beta_Q25) %>%
  arrange(-IQR) %>%
  .$Probe %>%
  .[1:10000]


mvals_variable_probes <-
  mvals_filtered[probenames_filtered %in% pca_probes_include, ] %>%
  apply(., 1, function(x) { x[is.na(x)] <- mean(x, na.rm = T); x }) %>%
  t

pca_results <-
  prcomp(mvals_variable_probes %>% t,
         center = T, scale = F)

pca_results <-
  prcomp(zvals_ruvg_cellcomp[ , probenames_filtered %in% pca_probes_include],
         center = T, scale = F)

pca_results_percentexp <-
  (pca_results$sdev^2) %>% `/`(., sum(.)) %>% `*`(100) %>% formatC(digits = 2) %>%
  gsub(" ", "", .) %>%
  paste0(" (", ., "%)") %>%
  paste0("PC", 1:length(.), .)

pca_results_df <-
  metadata_final_sample %>%
  bind_cols(pca_results$x[ , 1:5] %>%
              as_tibble %>% set_names(paste0("PC", 1:5)))



plot_pca_results <-
  function(varname, input_df = pca_results_df,
           input_percentexp = pca_results_percentexp,
           xaxis = "PC1", yaxis = "PC2") {
  if ((input_df[ , varname, drop = T] %>% is.numeric) |
      (input_df[ , varname, drop = T] %>% unique %>% length %>% `>`(., 6))) {
    ggplot(input_df,
           aes_string(xaxis, yaxis, color = varname)) +
      geom_point(alpha = 0.5) +
      xlab(pca_results_percentexp[1]) + ylab(pca_results_percentexp[2]) +
      theme_few()
  } else {
    ggplot(input_df,
           aes_string(xaxis, yaxis, color = varname, shape = varname)) +
      geom_point(alpha = 0.5) +
      xlab(input_percentexp[1]) + ylab(input_percentexp[2]) +
      theme_few()
  }
  }




# check PC assoc w/ various metadata -------------------------------------------
# & double check putative technical effects & nuisance effects (e.g., cell type)
# not associated with known features like sex/ancestry


plot_grid(
  plot_pca_results("CellTypePC1", xaxis = "PC1", yaxis = "PC2") +
  scale_color_viridis(name = "Cell Comp PC1", breaks = scales::pretty_breaks(n = 4)) +
    theme(legend.position = "bottom",
          legend.key.height = unit(0.2, "cm")),
  plot_pca_results("sex") +
    scale_color_manual(name = "Sex (male)", values = c("black", "blue")) +
    scale_shape_manual(name = "Sex (male)", values = c(1, 15)) + 
    theme(legend.position = "bottom", legend.key.height = unit(0.2, "cm")),
  plot_pca_results("genoPC1") +
    scale_color_viridis(name = "Genotype PC1", option = "D", breaks = scales::pretty_breaks(4)) +
    theme(legend.position = "bottom", legend.key.height = unit(0.2, "cm")),
  nrow = 1
)
  
ggsave("SupplFig2.png", width = 10, height = 3)
system('open SupplFig2.png')
  
plot_pca_results("EPIC_WellPairs") +
  facet_grid( ~  EPIC_WellPairs, labeller = as_labeller(function(x) { paste0("Row: ", x) })) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 0, lty = 3) +
  scale_shape(name = "EPIC Row") +
  scale_color_brewer(name = "EPIC Row", palette =  "Dark2") +
  theme(legend.position = "bottom")
  
ggsave("SupplFig2.png", width = 8, height = 2.75)
system('open SupplFig2.png')

plot_pca_results("CellTypePC1", xaxis = "PC1", yaxis = "PC2")
plot_pca_results("CellTypePC1", xaxis = "PC2", yaxis = "PC3")
plot_pca_results("EPIC_WellPairs", xaxis = "PC3", yaxis = "PC4")




plot_grid(
  plot_pca_results("CellTypePC1", xaxis = "PC1", yaxis = "PC2") +
    scale_color_viridis(name = "Cell Comp PC1", breaks = scales::pretty_breaks(n = 4)) +
    theme(legend.position = "bottom",
          legend.key.height = unit(0.2, "cm")),
  plot_pca_results("sex") +
    scale_color_manual(name = "Sex (male)", values = c("black", "blue")) +
    scale_shape_manual(name = "Sex (male)", values = c(1, 15)) + 
    theme(legend.position = "bottom", legend.key.height = unit(0.2, "cm")),
  plot_pca_results("genoPC1") +
    scale_color_viridis(name = "Genotype PC1", option = "D", breaks = scales::pretty_breaks(4)) +
    theme(legend.position = "bottom", legend.key.height = unit(0.2, "cm")),
  plot_pca_results("EPIC_WellPairs") +
    geom_hline(yintercept = 0, lty = 3) +
    geom_vline(xintercept = 0, lty = 3) +
    scale_shape(name = "EPIC Row") +
    scale_color_brewer(name = "EPIC Row", palette =  "Dark2") +
    theme(legend.position = "bottom"),
  nrow = 2
)

ggsave("SupplFig2.png", width = 8, height = 6)
system('open SupplFig2.png')




plot_pca_results("CellTypePC1") +
  scale_color_viridis_c() +
  facet_grid(~ Timepoint)
plot_pca_results("Timepoint") +
  facet_grid(~ Timepoint)

plot_pca_results("Timepoint") + facet_grid(~ sex)
plot_pca_results("Intervention") + facet_grid(~ Timepoint)
plot_pca_results("EPIC_WellPairs") + facet_grid(~ EPIC_Date_Scanned)
plot_pca_results("EPIC_WellPairs") + facet_grid(EPIC_WellPairs ~  .)
plot_pca_results("EPIC_Date_Scanned") +
  facet_grid(~ EPIC_Date_Scanned)
plot_pca_results("Specificity2")
plot_pca_results("EPIC_chip") 
plot_pca_results("sex")
plot_pca_results("Race")  + facet_grid(~ sex)
plot_pca_results("genoPC1") + scale_color_viridis_b()

plot_pca_results("outcome_wc_current", xaxis = "PC4", yaxis = "PC5") +
  scale_color_viridis_b() +
  facet_grid(~ Timepoint)
plot_pca_results("outcome_bmi_current", xaxis = "PC2", yaxis = "PC3") + scale_color_viridis_b() +
  facet_grid(~ Timepoint)
plot_pca_results("outcome_bmi_current", xaxis = "PC3", yaxis = "PC4") + scale_color_viridis_b() +
  facet_grid(~ Timepoint)
plot_pca_results("outcome_bmi_current", xaxis = "PC4", yaxis = "PC5") + scale_color_viridis_b() +
  facet_grid(~ Timepoint)

table(pca_results_df$Race, pca_results_df$EPIC_WellPairs)
table(pca_results_df$sex, pca_results_df$EPIC_WellPairs)

table(pca_results_df$EPIC_WellPairs, pca_results_df$sex)
table(pca_results_df$EPIC_WellPairs, pca_results_df$Race)

lm(PC1 ~ EPIC_row, pca_results_df) %>% anova()
lm(PC2 ~ EPIC_WellPairs, pca_results_df) %>% anova()

lm(PC1 ~ CellTypePC1, pca_results_df) %>% anova()
lm(PC1 ~ sex, pca_results_df) %>% anova()
lm(PC2 ~ genoPC1, pca_results_df) %>% anova()
lm(PC2 ~ race, pca_results_df) %>% anova()

lm(PC2 ~ outcome_bmi_current, pca_results_df) %>% anova()



dev.off()
