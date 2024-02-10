# ==============================================================================
# B/04_delta_correlations.R
# correlations between delta outcomes
# ==============================================================================





# load data, libraries ---------------------------------------------------------

load("Output/20220425/metadata_final.Rdata")

library(tidyverse)
library(tidygraph)
library(ggraph)
library(corrplot)



# check correlations between outcomes ------------------------------------------

delta_cor <-
  metadata_final_subject %>%
  select(contains("outcomedelta") & ( ends_with("30") | ends_with("60") )) %>%
  set_names(., names(.) %>% gsub("outcomedelta_", "Δ", .) %>% gsub("_", "", .)) %>%
  cor(use = "pairwise.complete.obs", method  = "spearman")

delta_cor[lower.tri(delta_cor, diag = T)] <- NA

delta_cor %>% .[ colnames(.) %>% grepl("Fat|weight|wc", .),
                 colnames(.) %>% grepl("Fat|weight|wc", .)]

# correlations between outcomes
edge_threshold <- 0.33 # 0.5
delta_cors <-
  as_tbl_graph(delta_cor, directed = F) %>%
  activate(edges) %>%
  filter(!is.na(weight)) %>%
  filter(abs(weight) > edge_threshold) %>%
  mutate(corlabel = if_else(abs(weight) > edge_threshold, formatC(weight, digits = 2), "")) %>%
  activate(nodes) 

# 800 x 600
set.seed(1234)
ggraph(delta_cors, layout = 'mds') +
  geom_edge_fan(aes(width = abs(weight), color = weight),#, label = corlabel),
                alpha = 0.8, label_size = 3) +
  geom_node_point(size = 3) +
  geom_node_text(aes(label = name), repel = T, size = 3, max.overlaps = 40) +
  scale_edge_colour_gradient2(name = "r", limits = c(-1, 1)) +
  theme_void(base_size = 10) +
  scale_edge_width(name = "|r|", range = c(0.25, 3))

