# ==============================================================================
# E/03_heatmap_caret_intervention.R
# plot multiple sites at a time
# e.g., rf/lasso selected sites for intervention effect
# ==============================================================================



load("Trt_DMethyl_RUVg_covarB.Rdata")

library(caret)
library(ComplexHeatmap)


lasso_variables <-
  lasso_fit$finalModel %>%
  varImp(scale = F) %>%
  filter(Overall != 0) %>%
  arrange(Overall) %>%
  row.names()

heatmap_to_plot <-
  caret_input_dat %>%
  select(all_of(lasso_variables))

set.seed(12345)
dcr_vs_imf_heatmap <-
  Heatmap(matrix = heatmap_to_plot %>% as.matrix,
        name = "3mo Delta\nMethylation",
        col = circlize::colorRamp2(c(-3, -2, -1, 0, 1, 2, 3), viridis::magma(7)),
        border = T,
        column_names_gp = gpar(fontsize = 6)#,
        ) +
  HeatmapAnnotation(Intervention = caret_input_dat$Intervention,
                col = list(Intervention = c(DCR = "#8c510a", IMF = "#5ab4ac")),
                border = T, show_annotation_name = c(F),
                annotation_width = 0.1,
                width = unit(0.1, "mm"),
                which = "row")

dcr_vs_imf_heatmap

