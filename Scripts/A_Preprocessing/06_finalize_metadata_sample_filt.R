# ==============================================================================
# A/06_finalize_metadata_sample_filt.R
# establish final sample-wise & subject-wise metadata
# filter by participant, omit outliers in cardiometabolic variables
# ==============================================================================




# metadata filtering -----------------------------------------------------------

filter_subjects_outliers <- T # * 
filter_subjects_dupes <- !grepl("_2", metadata$Sample_Name)
filter_subjects <- filter_subjects_outliers & filter_subjects_dupes
filter_samples_baselineonly <- metadata$Timepoint == "BL"

# * have considered omitting putative outliers:
# ID = "026", "039" (methyl), "055" (inferred cell comp), "003" (6mo outcomes)
# but retained due to small sample size, lack of independent justification

metadata_filtered <-
  metadata[filter_subjects, ]



# metadata filtered by sample (124 x 248) --------------------------------------

metadata_filtered_sample <-
  expand_grid(IID = metadata_filtered$IID %>% unique,
              Timepoint = c("BL", "3M")) %>%
  arrange(IID, rev(Timepoint)) %>%
  left_join(., metadata_filtered, by = c("IID", "Timepoint")) 



# filtered by subject (62 x 248) -----------------------------------------------

# note cell-type PCs vary by time, so BL & 3mo averaged for subject-level EDA;
# however, the sample-specific cell comp accounted for
# because of the residualization process
# (in script 03, EDA suggests far more variation between-patient than btwn-time)

mean_PCA_cellcomp <-
  metadata_filtered_sample %>% group_by(IID) %>%
  summarize(CellTypePC1 = mean(CellTypePC1), CellTypePC2 = mean(CellTypePC2))

metadata_filtered_subject <-
  metadata_filtered_sample %>%
  filter(!duplicated(IID)) %>% # grab the baseline values to de-dupe
  select(-CellTypePC1, -CellTypePC2) %>%
  left_join(., mean_PCA_cellcomp, by = "IID")




# remove extreme values from outcomes, outcome delta ---------------------------
# due to relatively small sample size, remove extreme values: elected to set NAs
# versus thresholding (e.g., <- quantile(x, 0.95)) bc of small sample size
# - 5sigma arbitrary threshold, implicitly implies bell curve like distribution 
#   (from EDA, most changes do look approx normal)
# - approach is likely conservative wrt to finding methylation associations
#   due to reduction in power from sample size, variance in outcome

omit_5sd_outlier <-
  function(y, multiplier = 5) {
    outmean <- mean(y, na.rm = T)
    outsd <- sd(y, na.rm = T)
    filter_is_outlier <-
      between(y, outmean - multiplier*outsd, outmean + multiplier*outsd) %>% `!`
    y[filter_is_outlier == T] <- NA
    y
  }

# final datasets created after omitting 5sd outliers
metadata_final_sample <-
  metadata_filtered_sample %>%
  mutate_at(
    names(.) %>% .[grepl("outcomedelta_|outcome_", .)],
    omit_5sd_outlier
  )
metadata_final_subject <-
  metadata_filtered_subject %>%
  mutate_at(
    names(.) %>% .[grepl("outcomedelta_|outcome_", .)],
    omit_5sd_outlier
  )

# check that final metadata differs (outlier-removed)
all_equal(metadata_final_sample, metadata_filtered_sample,
          ignore_col_order = F, ignore_row_order = F)
all_equal(t(metadata_final_sample), t(metadata_filtered_sample),
          ignore_col_order = F, ignore_row_order = F)




# final export -----------------------------------------------------------------

save(metadata_final_sample, metadata_final_subject,
     file = paste0("./Output/", date_export, "/metadata_final.Rdata"))

paste0("./Output/", date_export, "/metadata_final.Rdata")

