# ==============================================================================
# A/04_RUVg.R
# control probes from sesame (script 01) --> run RUVg
# ==============================================================================




library(tidyverse)
library(RUVSeq)

date_export <- "20220425"

ztransform <-
  function(x) { (x - mean(x)) / sd(x)}

# load control probe intensities -----------------------------------------------
# log2 & z-transform before input

# 1269 x 128
control_probes <-
  read_tsv("Output/20220421/methyl_controls.tsv", col_names = F) %>% 
  as.matrix() %>%
  `+`(., 1) %>% log2 %>% apply(., 1, ztransform)

ruvg_out <-
  RUVg(t(control_probes), k = 3, isLog = T)

ruvg_out$W %>%
  as_tibble() %>%
  bind_cols(Sample_Name = targets$Sample_Name, .) %>%
  write_tsv(file = paste0("./Output/", date_export, "/ruvg.tsv"))




# evaluation -------------------------------------------------------------------

# k = 3 is plateau of elbow plot

library(vegan)
map_dbl(1:10,
        function(k) {
          ruvg_out <-
            RUVg(t(control_probes), k = k, isLog = T)
          permanova_out <-
            adonis(control_probes ~ ruvg_out$W, method = "euclidean")
          permanova_out$aov.tab$R2[1]
        }) %>% plot


# further evaluation not shown: ------------------------------------------------

# * also iterative process with PCA in script "B/04_methyl_PCA.R",
# where looked for visual examination / associations of apparent batch effects
# (namely, diffuse/inconsistent effect direction by assay row) with first PCs
# * also evaluated COMBAT, RUVinv/RUVr: did not seem to remove apparent
# technical effect given same number of factors / conceptually less appropriate
# for the situation (not a unified row effect, multiple outcomes to residualize)

