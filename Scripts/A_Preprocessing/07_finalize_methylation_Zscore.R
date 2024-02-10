# ==============================================================================
# A/07_finalize_methylation_Zscore.R
# calculate Z-scores via residualization, standardization
# ==============================================================================




# load data --------------------------------------------------------------------

screen
qlogin -R rusage[mem=24]
cd $HOME/Analyses/DRIFT2_Methylation
module load R/4.1.0_beta
R

library(tidyverse)
library(data.table)
library(sesameData)

date_export <- "20220425"

# data.table of 865918 x 128 samples, col names = IID-Timepoint
methyl_betas <-
  fread("./Output/methyl_betas.tsv")

# 865918 Illumina probenames (e.g., "cg00025216")
probenames <-
  read_lines("./Output/20220421/probenames.txt")

# metadata from script 08
load("./Output/20220425/metadata_final.Rdata")



# load annotations -------------------------------------------------------------

# sesame betas sorted by probenames (cg000)
# the manifest sorted by chr -> position, sort by probename
# identical(probenames, EPIC.hg38.manifest_df$probe)

EPIC.hg38.manifest <- sesameDataGet('EPIC.hg38.manifest')
EPIC.hg38.manifest_df <-
  EPIC.hg38.manifest %>%
  as_tibble() %>%
  bind_cols(probe = EPIC.hg38.manifest@ranges@NAMES, .) %>%
  arrange(probe)



# criteria for probe inclusion -------------------------------------------------

filter_probes_autosomal <-
  EPIC.hg38.manifest_df$seqnames %in% paste0("chr", 1:22)
filter_probes_nonSNP <- 
  EPIC.hg38.manifest_df$probeType != "rs"

filter_probes_anyNA <-
  methyl_betas %>% apply(., 1, function(x) { is.na(x) %>% any} )
probes_beta_range <-
  apply(methyl_betas, 1, function(x) { max(x) - min(x) })
filter_probes_range_0.05 <-
  probes_beta_range > 0.05

# final probe filter (463,074 probes)
filter_probes <-
  !filter_probes_anyNA & filter_probes_range_0.05 &
  filter_probes_autosomal & filter_probes_nonSNP




# beta --> m-values ------------------------------------------------------------

methyl_betas_filtered <- copy(methyl_betas)
methyl_betas_filtered <- methyl_betas_filtered[filter_probes, ]

eps <- 1e-6

mvals_filtered <-
  copy(methyl_betas_filtered)
mvals_filtered <-
  mvals_filtered %>%
  .[, (names(.)) :=
      lapply(.SD, function(x) { log( (x + eps) / ( 1 - x + eps) ) } )]

dim(mvals_filtered) 
probenames_filtered <-
  probenames[filter_probes]

# should give data.table of 463963 probes x 124 samples 
# where samples are in order by IID, with baseline followed by 3mo value
mvals_filtered <- 
  mvals_filtered[ ,  metadata_final_sample$Sample_Name, with = F]

ztransform <-
  function(x) { (x - mean(x))/sd(x) }




# m-values  --> Z-scores -------------------------------------------------------
# yields standardized, covariate-adjusted m-values (Z-scores)

# adjusted for RUVg components
zvals_ruvg <-
  copy(mvals_filtered)
zvals_ruvg <-
  apply(zvals_ruvg, 1,
        function(y) { lm(y ~ RUV1 + RUV2 + RUV3,
                         metadata_final_sample) %>% resid %>% ztransform })

# adjusted for RUVg components & inferred cell-type PCs
zvals_ruvg_cellcomp <-
  copy(mvals_filtered)
zvals_ruvg_cellcomp <-
  apply(zvals_ruvg_cellcomp, 1,
        function(y) { lm(y ~ RUV1 + RUV2 + RUV3 + CellTypePC1 + CellTypePC2,
                         metadata_final_sample) %>% resid %>% ztransform })



# export -----------------------------------------------------------------------

# basic version, essentials only for modeling function
save(probenames_filtered, mvals_filtered, zvals_ruvg, zvals_ruvg_cellcomp,
     file = paste0("./Output/", date_export, "/final_methylation_modeling.Rdata"))

# some extended output for probewise summary stats, Z --> b conversion, etc
save(probenames_filtered, zvals_ruvg, zvals_ruvg_cellcomp,
     methyl_betas_filtered, mvals_filtered,
     file = paste0("./Output/", date_export, "/final_methylation_full.Rdata"))






