# ==============================================================================
# B/05_celltype_pcas.R
# check cell-type loadings, PCs
# ==============================================================================




# load data, libraries ---------------------------------------------------------

load("Output/20220425/metadata_final.Rdata")

library(tidyverse)
library(tidygraph)
library(ggraph)
library(corrplot)
library(ggbeeswarm)



# extract cell types, pca ------------------------------------------------------

# extract columns labeled "CellType_"
# estimated in pre-processing script A02
cellcounts_by_sample <-
  metadata_final_sample %>%
  select(IID, Timepoint, contains("CellType_"))

# check correlations by cell-type
# strong neu <--> lymphocyte inverse cor (seems captured in PC1),
# which in turn has moderate r = 0.3 - 0.6 to other celltypes (PC2)

cellcounts_by_sample %>%
  mutate(CellType_Lymphocytes =
           CellType_CD8T + CellType_CD4T + CellType_NK + CellType_Bcell) %>%
  select(contains("CellType_")) %>%
  cor(method = "spearman")

# run pca
pca_by_cellcounts <-
  cellcounts_by_sample %>%
  mutate_at(3:ncol(.), function(p) { log(p / (1 - p))}) %>%
  select(-IID, -Timepoint) %>%
  prcomp(scale = T)

pca_percents <-
  (pca_by_cellcounts$sdev ^ 2) %>%
  `/`(., sum(.)) %>%
  `*`(., 100) %>%
  formatC(digits = 3) %>%
  paste0("PC", 1:length(.), " (", ., "%)")

pca_loadings <-
  pca_by_cellcounts$rotation %>%
  .[ , 1:2] %>%
  as_tibble(rownames = "CellType") %>%
  mutate(PC1 = PC1 / max(PC1) *
           (range(metadata_final_sample$CellTypePC1) %>% abs %>% min),
         PC2 = PC2 / max(PC2) *
           (range(metadata_final_sample$CellTypePC2) %>% abs %>% min))




# show pca, colored by Neu -----------------------------------------------------

ggplot(metadata_final_sample,
       aes(CellTypePC1, CellTypePC2)) +
  geom_hline(lty = 3, alpha = 0.5, yintercept = 0) +
  geom_vline(lty = 3, alpha = 0.5, xintercept = 0) +
  geom_point(shape = 21, size = 2, aes(fill = CellType_Neu * 100)) +
  geom_segment(data = pca_loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.1, "inches"))) + # , color = "#b2182b"
  geom_text(data = pca_loadings,
            aes(x = if_else(PC1 < 0, PC1 - 0.3, PC1 + 0.3),
                y = PC2 + 0.3,
                label = gsub("CellType_", "", CellType)), size = 3) +
  # # potential outliers w.r.t. cell-type,
  # # but no clinical rationale for exclusion so do not omit
  # geom_text_repel(aes(label = if_else(CellTypePC2 > 3 | CellTypePC2 < -2.75, IID, "")),
  #                 size = 3,  color = "#b2182b") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  xlab(pca_percents[1]) + ylab(pca_percents[2]) +
  ggtitle("PCA on Estimated Cell %s") +
  scale_fill_viridis_c(name = "Neutrophil\nFreq (%)", direction = -1) +
  NULL

# 500 x 400
metadata_final_sample %>%
  select(contains("Celltype_")) %>%
  pivot_longer(cols = 1:ncol(.)) %>%
  transmute(CellType = gsub("CellType_", "", name),
            Percent = value * 100) %>%
  ggplot(data = ., aes(x = Percent, y = CellType)) +
  ylab("Cell Type") +
  scale_x_continuous(name = "Estimated Cell Frequency (%)", limits = c(0, 80)) +
  geom_density_ridges(jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
  theme_ridges()




# calculate delta cell-types ---------------------------------------------------
# visualizing if participants show shift after 3-mo w/in trial

deltacellcounts_by_sample <-
  metadata_final_sample %>%
  select(IID, Timepoint, contains("CellType_")) %>%
  pivot_longer(cols = 3:ncol(.), names_to = "CellType", values_to = "Percent") %>%
  pivot_wider(id_cols = IID, names_from = c("CellType", "Timepoint"), values_from = Percent) %>%
  transmute(IID,
            CD8T = CellType_CD8T_3M - CellType_CD8T_BL, CD4T = CellType_CD4T_3M - CellType_CD4T_BL,
            NK = CellType_NK_3M - CellType_NK_BL, Bcell = CellType_Bcell_3M - CellType_Bcell_BL,
            Mono = CellType_Mono_3M - CellType_Mono_BL, Neu = CellType_Neu_3M - CellType_Neu_BL) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "CellType", values_to = "Percent") %>%
  mutate(CellType = gsub("CellType_", "", CellType))

ggplot(data = deltacellcounts_by_sample,
       aes(x = Percent, y = CellType)) +
  scale_x_continuous(name = "3-Month Change in Cell Frequency (%)") +
  geom_density_ridges(jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) +
  theme_ridges() + ylab("Cell Type") +
  geom_vline(xintercept = 0, lty = 3)

# decreases in neutrophils
deltacellcounts_by_sample %>%
  group_by(CellType) %>%
  summarize(MeanDelta = mean(Percent))



# compare deltas by intervention -----------------------------------------------
# observe some mean differences, but extremely high variation &
# unclear DNAme measurement, intra-individual noise

deltacellcounts_by_sample %>%
    left_join(., metadata_final_subject, by = c("IID")) %>%
  group_by(CellType, Intervention) %>%
  summarize(Percent = summaryMeanSD(Percent, digits = 3)) %>%
  pivot_wider(names_from = Intervention, values_from = Percent)

deltacellcounts_by_sample %>%
  left_join(., metadata_final_subject, by = c("IID")) %>%
  ggplot(data = ., aes(Intervention, Percent)) +
  geom_quasirandom(alpha = 0.24) +
  facet_grid(~ CellType) +
  ggthemes::theme_few() + 
  stat_summary() +
  geom_hline(yintercept = 0, lty = 3)




