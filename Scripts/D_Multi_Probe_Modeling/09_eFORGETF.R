# ==============================================================================
# D/09_eFORGETF.R
# p-value correct eFORGE-TF analyses
# ==============================================================================



library(tidyverse)

eforge_res <-
  read.table("./Output/20220425/eFORGE-TF/eFORGE-TF.2022-04-27T08_52_56.822Z-Trt_DMethyl_setA.csv", skip = 7, encoding = "UTF-16LE", skipNul = T) %>%
  arrange(V3)
eforge_res$FDR <-
  p.adjust(eforge_res$V3, "fdr")

eforge_res


