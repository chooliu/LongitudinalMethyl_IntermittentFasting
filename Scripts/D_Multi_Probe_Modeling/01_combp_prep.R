# ==============================================================================
# D/01_combp_prep.R
# format results tables from modeling output (scripts in C) into .bed
# ==============================================================================




# below run interactively ------------------------------------------------------

screen -S combp
qlogin -R rusage[mem=8]


mkdir $HOME/Analyses/DRIFT2_Methylation/Output/20220425/combp
cd $HOME/Analyses/DRIFT2_Methylation/

module load R/4.1.0_beta
R



# fread results --> change names --> export as .bed ----------------------------

library(data.table)
library(sesameData)
library(tidyverse)

# sesame annotations for probe <--> hg38 location
EPIC.hg38.manifest <-
  sesameDataGet('EPIC.hg38.manifest')
EPIC.hg38.manifest

bedoutputdir <- "./Output/20220425/combp/bedfiles/"
dir.create("./Output/20220425/combp/")
dir.create(bedoutputdir)

# target results .tsv from scripts in part C
filepaths <-
  list.files(path = "./Output/20220425/ResultsTables", full.names = T)

# e.g., if only need to run subset of results tables
# e.g., for the baseline --> short term change in outcomes
# filepaths <- filepaths %>% .[grepl("Methyl0", .)] %>%  .[grepl("_30|_60", .)]

# run filepaths
sapply(filepaths,
       function(filepath_results_tsv) {

         results <- fread(filepath_results_tsv)
         results <- results[!is.na(bacon_Pz), ]
         results_ranges <-
           EPIC.hg38.manifest[results$Probe] %>%
           as_tibble %>%
           bind_cols(p = results$bacon_Pz) %>%
           transmute(`#chrom` = seqnames %>% as.character, start, end, p) %>%
           arrange(`#chrom`, start)
         
         write_tsv(results_ranges,
                   file = paste0(bedoutputdir,
                                 gsub("tsv", "bed", basename(filepath_results_tsv))))
         
         paste0(bedoutputdir, gsub("tsv", "bed", basename(filepath_results_tsv)))
       }
)

