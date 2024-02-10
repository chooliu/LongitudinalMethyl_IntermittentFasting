# ==============================================================================
# C/04_annotate_DMPs.R 
# compile information about each probe to aid in interpretation of DMP results
# ==============================================================================



# load data --------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(sesameData)

load("Output/20220425/final_methylation_full.Rdata")
date_export <- "20220425"



# load EPIC annotations --------------------------------------------------------

# sesame betas sorted by probenames (cg000)
# the manifest sorted by chr -> position, sort by probename
# identical(probenames, EPIC.hg38.manifest_df$probe)

EPIC.hg38.manifest <- sesameDataGet('EPIC.hg38.manifest')
EPIC.hg38.manifest_df <-
  EPIC.hg38.manifest %>%
  as_tibble() %>%
  bind_cols(probe = EPIC.hg38.manifest@ranges@NAMES, .) %>%
  arrange(probe)

# https://zwdzwd.github.io/InfiniumAnnotation
EPIC.hg38.manifest_additional <-
  fread("./Output/20220425/EPIC.hg38.manifest.gencode.v36.tsv") %>%
  transmute(probe = probeID, DistToTSS = distToTSS, CGType = CGIposition, ENST = transcriptIDs) %>%
  as_tibble()

EPIC.hg38.manifest_df <-
  EPIC.hg38.manifest_df %>%
  left_join(EPIC.hg38.manifest_additional, by = "probe")



# summary stats on betas -------------------------------------------------------

probes_beta_quantiles <-
  apply(methyl_betas_filtered, 1, function(x) { quantile(x, na.rm = T) })

probe_summarystats <-
  probes_beta_quantiles %>%
  t() %>%
  as_tibble() %>%
  set_names(c("Beta_Min", "Beta_Q25", "Beta_Median", "Beta_Q75", "Beta_Max")) %>%
  mutate(Beta_Range = Beta_Max - Beta_Min) %>%
  bind_cols(Probe = probenames_filtered)

probe_EPIC_annot <-
  EPIC.hg38.manifest_df %>%
  filter(probe %in% probenames_filtered) %>%
  transmute(Probe = probe,
            Chr = seqnames,
            Start = start, End = end, Strand = strand,
            Gene = gene_HGNC, ENST,
            DistToTSS, CGType)
levels(probe_EPIC_annot$Chr) <- 
  probe_EPIC_annot$Chr %>%
  levels %>%
  gsub("chr", "", .)



# Z-value +/- 1 <--> estimated beta conversion ---------------------------------

conversion_factor <-
  sapply(1:nrow(mvals_filtered),
         function(i) {
           Fz <- ecdf(zvals_ruvg[, i])
           quantiles_Z <- Fz(c(-1, 0, 1))
           bvec <- methyl_betas_filtered[i, ] %>% as.numeric
           bl <- quantile(bvec, quantiles_Z) %>% as.numeric
           mean(c(bl[3] - bl[2], bl[2] - bl[1]))
         })

conversion_factor_cellcomp <-
  sapply(1:nrow(mvals_filtered),
         function(i) {
           Fz <- ecdf(zvals_ruvg_cellcomp[, i])
           quantiles_Z <- Fz(c(-1, 0, 1))
           bvec <- methyl_betas_filtered[i, ] %>% as.numeric
           bl <- quantile(bvec, quantiles_Z) %>% as.numeric
           mean(c(bl[3] - bl[2], bl[2] - bl[1]))
         })

# separate conversion factor for each z-score
# (+/- cell composition changes)
# but turned out to be very similar
cor(conversion_factor, conversion_factor_cellcomp) # 0.982
cor(conversion_factor, conversion_factor_cellcomp, method = 'spearman') # 0.993



# compile info, export ---------------------------------------------------------

probe_annotations <-
  left_join(probe_EPIC_annot, probe_summarystats, by = "Probe") %>%
  bind_cols(BetaConversion1 = conversion_factor,
            BetaConversion2 = conversion_factor_cellcomp)

save(probe_annotations,
     file = paste0("./Output/", date_export, "/probe_annotations.Rdata"))


