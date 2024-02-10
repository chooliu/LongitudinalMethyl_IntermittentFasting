# ==============================================================================
# B/02_baseline_correlations.R
# check baseline correlations between cell comp, ruv, demographics, outcomes
# ==============================================================================




# load data, libraries ---------------------------------------------------------

load("Output/20220425/metadata_final.Rdata")

library(tidyverse)
library(tidygraph)
library(ggraph)
library(corrplot)



# calculate correlations, set heatmap order ------------------------------------

# subset columns of interest
corplot_input <-
  metadata_final_subject %>%
  select(contains("CellType"), contains("_0"), contains("RUV"),
         sex, contains("geno"), age, Intervention) %>%
  mutate(Intervention = as.numeric(Intervention == "IMF"))

corplot_var_order <-
  c(names(corplot_input) %>% .[!grepl("outcome", .) & !grepl("CellType", .)] %>% sort,
    names(corplot_input) %>% .[grepl("CellType", .)] %>% sort,
    names(corplot_input) %>% .[grepl("outcome", .) & !grepl("CellType", .)] %>% sort)

corplot_input <-
  corplot_input[ , corplot_var_order] %>%
  mutate_all(as.double) %>%
  as.matrix

colnames(corplot_input) <-
  colnames(corplot_input) %>%
  tolower %>% gsub("outcome_|_0|celltype_", "", .)
row.names(corplot_input) <- names(corplot_input) 

# spearman, with asterisk on tiles w abs(cor) > 0.33
cor_mat <- cor(corplot_input, method = 'spearman', use = "complete") 
cor_asterisk <- abs(cor_mat) <= 0.33



# plot heatmap -----------------------------------------------------------------

png("FigSX_heatmap.png", width = 3000, height = 3000, res = 300)
  corrplot(cor_mat,
           p.mat = cor_asterisk,
           is.corr = T,
           tl.col = "black", 
           type = 'lower',
           sig.level = 0.05,
           method = 'square',
           insig = "label_sig",
           pch.cex = 0.8,
           diag  = T) %>%
  corrRect(name = c('age', 'bcell', 'bmi', 'weight'))
dev.off()



# alternative network viz, correlations btwn baseline values -------------------

baseline_cor <-
  metadata_final_subject %>%
  select(ends_with("0") & !contains("outcomedelta")) %>%
  set_names(., names(.) %>% gsub("outcome|_", "", .)) %>%
  cor(use = "pairwise.complete.obs", method  = "spearman")

baseline_cor[lower.tri(baseline_cor, diag = T)] <- NA

baseline_cor <-
  as_tbl_graph(baseline_cor, directed = F) %>%
  activate(edges) %>%
  filter(!is.na(weight)) %>%
  filter(abs(weight) > 0.33) %>%
  mutate(corlabel = if_else(abs(weight) > 0.33, formatC(weight, digits = 2), "")) %>%
  activate(nodes)

# 800 x 600
ggraph(baseline_cor, layout = 'auto') +
  geom_edge_fan(aes(width = abs(weight), color = weight),
                alpha = 0.8, label_size = 3) +
  geom_node_point(size = 3) +
  geom_node_text(aes(label = name), repel = T, size = 3, max.overlaps = 15) +
  scale_edge_colour_gradient2(name = "r", limits = c(-1, 1)) +
  theme_void() +
  scale_edge_width(name = "|r|",
                   range = c(0.5, 1.5)) +
  ggtitle("Cross-Sectional (Baseline) Correlations")


