# ==============================================================================
# A/02_estimate_cellcomp.R
# .idat --> beta values
# ==============================================================================




# run .idat targets (targets via script 01), -----------------------------------
# but using minfi infrastructure for compatibility

library(tidyverse)
library(minfi)

# (ideally would have one set of betas used across all analyses, but 
# since cell deconv reference dataset used minfi, this prob more accurate &
# the cell proportions are not of interest in and of themselves)

# load idats via minfi
# note: targets file in  load("./Output/20220421/sesame_output.Rdata")
rgset_drift <-
  read.metharray.exp(
    targets = targets %>%
      transmute(Basename = file_path_prefix))



# experimental hub for whole blood ref dataset (measured via EPIC) -------------

library(ExperimentHub)
hub <- ExperimentHub()  
query(hub, "FlowSorted.Blood.EPIC")  
FlowSorted.Blood.EPIC <- hub[["EH1136"]]  




# run & save cell counts -------------------------------------------------------

cellcomp <-
  estimateCellCounts2(rgset_drift,
                      compositeCellType = "Blood",   
                      processMethod = "preprocessNoob",  
                      cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),  
                      referencePlatform = "IlluminaHumanMethylationEPIC",  
                      IDOLOptimizedCpGs = IDOLOptimizedCpGs,
                      lessThanOne = T)

# export
cellcomp$counts %>%
  as_tibble(rownames = "Plate_ID") %>%
  write_tsv(file = paste0("./Output/", date_export, "/cellcomp.tsv"))
