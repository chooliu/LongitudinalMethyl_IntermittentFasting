# ==============================================================================
# B/03_check_genoPCA.R
# check genotyping PCA, compare to participant self-report 
# ==============================================================================




# load data, libraries ---------------------------------------------------------

load("Output/20220425/metadata_final.Rdata")

library(tidyverse)
library(RColorBrewer)



# load eigenvalues -------------------------------------------------------------

geno_eigenvals <-
  read_lines("./data/20201007-eigenstrat/ADA_eigenvals.out") %>%
  as.numeric()

genoPCA_percent_exp <-
  geno_eigenvals %>% `/`(., sum(.)) %>% .[1:10] 

# elbow plot, 2-3 genotype PCs
ggplot(data = NULL,
       aes(x = 1:10, y = genoPCA_percent_exp * 100)) +
  geom_point() + geom_line() +
  theme_few() +
  scale_x_continuous(breaks = 1:10) +
  xlab("PC") +
  ylab("Genotype PCA,\n% Variance Explained")
  
genoPCA_percent_exp <-
  genoPCA_percent_exp %>%
  `*`(., 100) %>% formatC(digits = 2) %>% paste0(" (", ., "%)") %>%
  paste0("Genotype PC", 1:10, .)




# plot compared to SIRE --------------------------------------------------------

# generally segregates with participant self-reports (SIRE)
# SIRE largely consistent with genoPCs

# PC1: black/african-american versus other ancestry
# PC2: hispanic versus non-Hispanic
# PC3: variability among black/african-american participants

palette_SIRE_color <-
  brewer.pal(n = 5, name = "Set1")

plot_grid(
  ggplot(metadata_final_subject,
         aes(genoPC1, genoPC2,
             shape = Ethnicity, color = Race)) +
    geom_point(size = 2, alpha = 0.5) +
    theme_few() +
    xlab(genoPCA_percent_exp[1]) + ylab(genoPCA_percent_exp[2]) +
    scale_color_manual(values = palette_SIRE_color) +
    scale_shape_manual(values = c(0, 19)) +
    theme(legend.position = "none"),
  
  ggplot(metadata_final_subject,
         aes(genoPC1, genoPC3,
             shape = Ethnicity, color = Race)) +
    geom_point(size = 2, alpha = 0.5) +
    theme_few() +
    xlab(genoPCA_percent_exp[1]) + ylab(genoPCA_percent_exp[3]) +
    scale_color_manual(values = palette_SIRE_color) +
    scale_shape_manual(values = c(0, 19)) +
    theme(legend.spacing.y = unit(0.1, "cm")),
  rel_widths = c(4, 6))


