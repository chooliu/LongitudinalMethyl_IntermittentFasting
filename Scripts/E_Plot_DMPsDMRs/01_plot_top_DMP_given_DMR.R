# ==============================================================================
# E/01_plot_top_DMP_given_DMR.R
# given DMR of interest, select & plot top effect size DMP from DMR
# either z-scores or raw beta
# ==============================================================================



# prep -------------------------------------------------------------------------

library(readxl)
library(ggthemes)
library(data.table)


# # raw data
# betas_raw <- fread("./Data/20220421-methylation/methyl_betas.tsv")
# probenamesall <- read_lines("./Data/20220421-methylation/probenames.txt")

# # corrected z-scores
# load("ModelInput/final_methylation_full.Rdata")
# load("ModelInput/20210512_metadata_final.Rdata")


load("./Output/20220425/metadata_final.Rdata")
( dmrs <- read_excel("TargetsToPlot/v2_DRIFT_topDMRs_forgraphs_1-30-23.xlsx", sheet = 1) )


library(sesameData)
EPIC.hg38.manifest <- sesameDataGet('EPIC.hg38.manifest')
EPIC.hg38.manifest_df <-
  EPIC.hg38.manifest %>%
  as_tibble() %>%
  bind_cols(probe = EPIC.hg38.manifest@ranges@NAMES, .) %>%
  arrange(probe)

# make concordant with outcome names given by SB
metadata_final_sample$outcomedelta_FatPercWB_60 <-
  metadata_final_sample$outcome_FatPercWB_6 - metadata_final_sample$outcome_FatPercWB_0

metadata_final_sample$outcomedelta_LDL_60 <-
  metadata_final_sample$outcomedelta_ldl_60




# helper fxn for plots ---------------------------------------------------------

make_dmp_bl_scatterplot <- 
  function(probe_in, outcome_in, xlab = "", ylab = "") {
    metadata_final_sample %>%
      bind_cols(probe = zvals_ruvg_cellcomp[ , probenames_filtered == probe_in]) %>%
      # # or raw beta data
      # bind_cols(probe = methyl_betas_filtered[ probenames_filtered == probe_in , , drop = T] %>%
      #             as_tibble %>% .[ , metadata_final_sample$Sample_Name] %>% unlist) %>%
      filter(Timepoint == "BL") %>%
      ggplot(data = ., aes_string(y = outcome_in, x = "probe")) +
      geom_point(alpha = 0.4, color = "#1d91c0") +
      theme_few() +
      scale_y_continuous(name = ylab) +
      scale_x_continuous(name = xlab) 
  }

( dmps <- read_excel("TargetsToPlot/v2_DRIFT_topDMRs_forgraphs_1-30-23.xlsx", sheet = 1) )
dmps <-
  dmps %>%
  mutate(probe_in = `Top DMP`,
         outcome_in = `Clinical Outcome` %>% str_sub(2, 1000) %>% paste0("outcomedelta_", ., "_", `Time Point`, "0"),
         #ylab = paste0(),
         xlab = paste0(`Top DMP`, " (β), ", str_split(`DMR Annotation`, pattern = "\\(", simplify = T) %>% .[1]) )
  
dmps <- 
  left_join(dmps,
            tibble(outcome_in = dmps$outcome_in %>% unique,
                   ylab = c("3mo Δ weight (lb)", "6mo Δ weight (lb)",
                            "6mo kg fat/kg whole body (%)", "3mo Δ triglycerides (mg/dl)",
                            "6mo Δ triglycerides (mg/dl)", "3mo Δ waist circumference (cm)",
                            "6mo Δ waist circumference (cm)", "6mo Δ diastolic (mm Hg)",
                            "6mo Δ systolic (mm Hg)", "6mo Δ insulin (uIU/ml)",
                            "6mo Δ LDL (uIU/ml)", "6mo Δ BMI (kg/m2)",
                            "3mo Δ glucose (mg/dl)", "6mo Δ leptin (ng/nl)",
                            "6mo Δ CRprotein (mg/l)", "6mo Δ LDL (uIU/ml)")),
            by = "outcome_in"
  )


dmp_plots <- 
  dmps %>%
  select(probe_in, outcome_in, xlab, ylab) %>%
  pmap(make_dmp_bl_scatterplot)


map(1:length(dmp_plots),
    ~ ggsave(filename = paste0("DMRtopprobe_zvals_", .x, ".png"),
             plot = dmp_plots[[.x]], width = 3, height = 2.5, dpi = 300))

map(1:length(dmp_plots),
    ~ ggsave(filename = paste0("DMRtopprobe_betas_", .x, ".png"),
             plot = dmp_plots[[.x]], width = 3, height = 2.5, dpi = 300))

         












