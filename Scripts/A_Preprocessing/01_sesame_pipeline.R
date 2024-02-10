# ==============================================================================
# A/01_sesame_pipeline.R
# .idat --> beta values
# ==============================================================================




# server setup -----------------------------------------------------------------
# (cmd line)

export PATH=/beevol/home/lincuini/Software/:$PATH
export LD_LIBRARY_PATH=/beevol/home/lincuini/Software/libpng-1.6.37/lib/:$LD_LIBRARY_PATH


screen
qlogin -R rusage[mem=48]
cd $HOME/Analyses/DRIFT2_Methylation
module load R/4.1.0_beta
R



# libraries --------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(readxl)
library(data.table)
library(Cairo)
library(sesame)
library(khroma)

# sample metadata --------------------------------------------------------------

date_export <- "20220421"
dir.create("Output")
dir.create(paste0("Output/", date_export))

folder_raw_dat <-
  paste0(Sys.getenv("HOME"), "/Data/DRIFT2/Methylation/")

samplesheet <-
  left_join(
    read_csv(
      paste0(folder_raw_dat, "SarahBorengasser_128_samplesheet_09052019new.csv"),
      skip = 7) %>%
      select(Sample_Name, Sample_Well, Sample_Plate),
    read_excel(
      paste0(folder_raw_dat, "Sarah_Borengasser_09042019Controls.xlsx")),
    by = c("Sample_Name" = "Sample Name")
    )

# manual entry from O.C. (Genomics Shared Resource) e-mail
# Thursday, October 8, 2020 3:39 PM
assay_scan_date <-
  tibble(EPIC_Date_Scanned =
           c(rep("09042019", 6), rep("09052019", 4), rep("09062019", 6)),
         EPIC_chip = c("203784950100", "203784950110", "203806760084", "203806760086", "203806760118", "203806760124",
                        "203784950152", "203806760090", "203806760119", "203806760153",
                        "203784950101", "203784950111", "203784950132", "203784950146", "203784950151", "203806760151"))

# file path example format:
# Slide_1 _203784950100/203784950100_R01C01_Grn.idat
# use to extract assay position (slide x row) information
targets <-
  tibble(file_path =
           list.files(path = folder_raw_dat, pattern = "*.idat|*.IDAT",
                      full.names = F, recursive = T)) %>%
  mutate(EPIC_chip = str_split(file_path, pattern = "/|_") %>% map_chr( ~ .[4]),
         EPIC_row = str_split(file_path, pattern = "/|_") %>% map_chr( ~ .[5])) %>%
  left_join(samplesheet,
            by = c("EPIC_chip" = "Sentrix Barcode",
                   "EPIC_row" = "Sentrix Position")) %>%
  mutate(file_path_prefix = gsub("_Red.idat|_Grn.idat", "", file_path) %>%
           paste0(folder_raw_dat, .)) %>%
  filter(!duplicated(file_path_prefix)) %>% # duplicated prefix due to Red/Grn
  left_join(assay_scan_date, by = "EPIC_chip") %>%
  mutate(EPIC_WellPairs = EPIC_row %>% as.factor() %>% # technical effect by row
           fct_collapse(`1/2` = c("R01C01", "R02C01"),
                        `3/4` = c("R03C01", "R04C01"),
                        `5/6` = c("R05C01", "R06C01"),
                        `7/8` = c("R07C01", "R08C01")))




# load .idat values ------------------------------------------------------------

sesameDataCacheAll()

methylation_raw <- 
  map(targets$file_path_prefix, readIDATpair)

methylation_betas <-
  targets$file_path_prefix %>%
  openSesame(pval.threshold	= 0.05)

probenames <- rownames(methylation_betas)

methylation_betas <-
  methylation_betas %>% 
  as.data.table %>%
  set_names(., targets$Sample_Name)

# quality control
inferred_sex <-
  map_chr(methylation_raw, inferSex)
inferred_sexkarotype <-
  map_chr(methylation_raw, inferSexKaryotypes)
inferred_ethnicity <-
  map_chr(methylation_raw, inferEthnicity)
inferred_age <-
  map_dbl(methylation_betas, ~ unlist(.x) %>%
            set_names(probenames) %>% predictAgeHorvath353())

sesame_qc <-
  methylation_raw %>%
  map(sesameQC) %>%
  map_dfr( ~ as_tibble(.x))

sesame_qc <-
  sesame_qc %>%
  bind_cols(InferredAge = inferred_age,
            InferredEthnicity = inferred_ethnicity,
            InferredSex = inferred_sex,
            InferredKarotype = inferred_sexkarotype)





# export resulting methylation values ------------------------------------------

write_lines(probenames,
            file = paste0("./Output/", date_export, "/probenames.txt"))

fwrite(
  methylation_betas,
  file = paste0("./Output/", date_export, "/methyl_betas.tsv"),
  sep = "\t")

write_tsv(
  sesame_qc,
  file = paste0("./Output/", date_export, "/sesame_qc.tsv"))

write_tsv(
  targets,
  file = paste0("./Output/", date_export, "/targets.tsv"))

controls_G <-
  methylation_raw %>% map_dfc(~ controls(.x)[, "G"])
controls_R <-
  methylation_raw %>% map_dfc(~ controls(.x)[, "R"])
write_tsv(
  bind_rows(controls_G,
        controls_R),
  file = paste0("./Output/", date_export, "/methyl_controls.tsv"),
  col_names = F
)

methylation_raw %>% map_dfc(~ .x@ctl[ , "R"]) %>% write_tsv(
  path = paste0("./Output/", date_export, "/controls_R.tsv")
)

save(
  methylation_raw, methylation_betas, targets, sesame_qc,
  file = paste0("./Output/", date_export, "/sesame_output.Rdata"))

  


# cursory plotting to check batch effects --------------------------------------
# 8 samples per plate, 16 plates

palette_colors <-
  c(colour("light")(8), colour("muted")(8))

export_expression_boxplots <-
  function(dat, fileout, ylabel) {
    par(mar = c(1, 1, 1, 1))
    CairoPNG(filename = paste0("./Output/", date_export, "/", fileout),
             dpi = 300, width = 1600, height = 1600)
    range <- dat %>% range(na.rm = T)
    layout(matrix(c(1, 1, 1,
                    2, 2, 2), nrow = 2, byrow = TRUE))
    boxplot(dat[ , 1:(8*8)],
            boxfill = rep(palette_colors[1:8], each = 8),
            xaxt = "n", ylim = range, ylab = ylabel)
    boxplot(dat[ , (8*8 + 1):(8*16)],
            boxfill = rep(palette_colors[9:16], each = 8),
            xaxt = "n", ylim = range, ylab = ylabel)
    dev.off()
  }

export_expression_boxplots(methylation_raw %>% map_dfc(~ .x$UG),
                           "rawintensities_UG.png", "UG")
export_expression_boxplots(methylation_raw %>% map_dfc(~ .x$MR),
                           "rawintensities_MR.png", "MR")
export_expression_boxplots(methylation_raw %>% map_dfc(~ .x$UG),
                           "rawintensities_UG.png", "UG")
export_expression_boxplots(methylation_raw %>% map_dfc(~ .x$UR),
                           "rawintensities_UR.png", "UR")




# calculate m-values -----------------------------------------------------------

eps <- 1e-6
methylation_mvals <- copy(methylation_betas)
methylation_mvals <-
  methylation_mvals %>%
  .[, (names(.)) :=
      lapply(.SD, function(x) { log( (x + eps) / ( 1 - x + eps) ) } )] %>%
  as.matrix

export_expression_boxplots(
  methylation_betas %>% as.matrix,
  paste0("./Output/", date_export, "/processed_betas.png"),
  "SeSAMe-Normalized Betas")
export_expression_boxplots(
  methylation_mvals,
  paste0("./Output/", date_export, "/processed_mvals.png"),
  "SeSAMe-Normalized M-Values")

CairoPNG(filename = paste0("./Output/", date_export, "/example_densities_betas.png"),
         dpi = 300, width = 1000, height = 1000)
par(cex = 0.5)
minfi::densityPlot(methylation_betas[ , 1:48] %>% as.matrix,
            sampGroups = targets$EPIC_chip[1:48],
            pal = palette_colors[1:5])
dev.off()

CairoPNG(filename = paste0("./Output/", date_export, "/example_densities_mvals.png"),
         dpi = 300, width = 1000, height = 1000)
par(cex = 0.5)
minfi::densityPlot(methylation_mvals[ , 1:48] %>% as.matrix,
            sampGroups = targets$EPIC_chip[1:48],
            pal = palette_colors[1:5])
dev.off()

