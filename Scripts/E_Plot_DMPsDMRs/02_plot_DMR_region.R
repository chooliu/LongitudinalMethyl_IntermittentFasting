# ==============================================================================
# E/02_plot_DMR_region.R
# given DMR of interest, plot all probes in window
# ==============================================================================



# prep -------------------------------------------------------------------------

library(readxl)
library(ggthemes)

# betas_raw <- fread("./Data/20220421-methylation/methyl_betas.tsv")
# probenamesall <- read_lines("./Data/20220421-methylation/probenames.txt")

# load("ModelInput/final_methylation_full.Rdata")

# load("ModelInput/20210516_probe_annotations.Rdata")
load("ModelInput/20210512_metadata_final.Rdata")
( dmrs <- read_excel("~/Downloads/v2_DRIFT_topDMRs_forgraphs_1-30-23.xlsx", sheet = 1) )


make_bl_dmrplot <- 
  function(dmr_chr, dmr_left, dmr_right, outcome_in, padding,
           custom_outcome_name = "Outcome",
           xlab = paste0(dmr_chr, ":", dmr_left, "-", dmr_right, "\n(length = ", dmr_right - dmr_left, ")" ),
           ylab = "Fraction Methylated (β)", nbreaks = 5) {
    
    probe_in_df <-
      EPIC.hg38.manifest_df %>% 
      filter(seqnames == dmr_chr, between(start, dmr_left - padding, dmr_right  + padding)) %>%
      # filter(probe %in% probenames_filtered) %>%
      transmute(probe, pos = start, indmr =  between(pos, dmr_left, dmr_right))
    
    # may be less visually confusing to discretize
    categ_outcome <-
      metadata_final_sample[ , outcome_in, drop = T] # %>%
      # cut_interval(., n = nbreaks) %>%
      # `levels<-`(c(1:nbreaks))
    
    metadata_final_sample$categ_outcome <- categ_outcome
    
    probe_val <-
      # zvals_ruvg_cellcomp[ , probenames_filtered %in% probe_in_df$probe]
      # methyl_betas_filtered[ probenames_filtered %in% probe_in_df$probe, ] %>%
      #             as_tsibble %>% .[ , metadata_final_sample$Sample_Name] %>%
      # t
      betas_raw[ probenamesall %in% probe_in_df$probe, ] %>%
      as_tibble %>% .[ , metadata_final_sample$Sample_Name] %>%
      t
    
    tmp_plot <-
      metadata_final_sample %>%
      select(IID, categ_outcome, Timepoint) %>%
      bind_cols(., probe_val %>% as_tibble %>% set_names(probe_in_df$probe)) %>%
      # filter(categ_outcome %in% c(1, 4)) %>%
      # group_by(IID, categ_outcome, Timepoint) %>%
      group_by(categ_outcome, Timepoint) %>%
      summarize_if(is.numeric, function(x) { median(x, na.rm = T) }) %>%
      pivot_longer(cols = names(.) %>% .[grepl("cg", .)], names_to = "probe") %>%
      left_join(probe_in_df, by = "probe")
    
    tmp_plot %>%
      filter(Timepoint == "BL" & !is.na(categ_outcome)) %>%
      # ggplot(data = ., aes(x = pos, y = value,
      #                      group = IID, color = categ_outcome)) +
      ggplot(data = ., aes(x = pos, y = value,
                           group = categ_outcome, color = categ_outcome)) +
      # geom_rect(data = NULL, aes(xmin = dmr_left, xmax = dmr_right, ymin = -Inf, ymax = Inf),
      #           fill = "#ffebff", alpha = 0.10, color = "#ffebff") +
      # geom_point(alpha = 0.4) +
      # geom_smooth(se = F) +\
      geom_line() +
      theme_few() +
      # scale_color_manual(name = custom_outcome_name,
      #                    values = c(`1` = "#000000", `2` = "#6e016b", `3` = "#8c6bb1", `4` = "#8c96c6", `5` = "#9ebcda"),
      #                    labels = c("1 (high)", "2", "3", "4", "5 (low)")) +
      scale_color_binned(name = custom_outcome_name) + #,
                         # breaks = probe_in_df$pos, labels = if_else(probe_in_df$probe %in% probenames_filtered, "*", ""),
                         #limits = c(dmr_left - padding, dmr_right + padding )) +
      theme(legend.position = "bottom") +
      scale_y_continuous(name = ylab) 
      
  }



# plot few examples ------------------------------------------------------------

index <- 1
dmrs[index, ]
make_bl_dmrplot(dmrs$`DMR Chromosome`[index] %>% paste0("chr", .),
                dmrs$`DMR Start`[index] %>% as.numeric(),
                dmrs$`DMR End`[index] %>% as.numeric(),
                dmrs$`Clinical Outcome`[index] %>% str_sub(., 2, 1000) %>%
                  #gsub(., "FatPercWB", "WBfatperc") %>%
                  paste0("outcomedelta_", .,  "_", dmrs$`Time Point`[index], "0"),
                "3mo ΔWeight (kg)", padding = 100, nbreaks = 5
)


index <- 2
dmrs[index, ]
make_bl_dmrplot(dmrs$`DMR Chromosome`[index] %>% paste0("chr", .),
                dmrs$`DMR Start`[index] %>% as.numeric(),
                dmrs$`DMR End`[index] %>% as.numeric(),
                dmrs$`Clinical Outcome`[index] %>% str_sub(., 2, 1000) %>%
                  #gsub(., "FatPercWB", "WBfatperc") %>%
                    paste0("outcomedelta_", .,  "_", dmrs$`Time Point`[index], "0"),
                "Weight (kg)", padding = 100, nbreaks = 5
)


index <- 3
dmrs[index, ]
make_bl_dmrplot(dmrs$`DMR Chromosome`[index] %>% paste0("chr", .),
                dmrs$`DMR Start`[index] %>% as.numeric(),
                dmrs$`DMR End`[index] %>% as.numeric(),
                dmrs$`Clinical Outcome`[index] %>% str_sub(., 2, 1000) %>%
                  gsub(., "FatPercWB", "WBfatperc") %>%
                  paste0("outcomedelta_", .,  "_", dmrs$`Time Point`[index], "0"),
                "Weight (kg)", padding = 100, nbreaks = 5
)

index <- 4
dmrs[index, ]
make_bl_dmrplot(dmrs$`DMR Chromosome`[index] %>% paste0("chr", .),
                dmrs$`DMR Start`[index] %>% as.numeric(),
                dmrs$`DMR End`[index] %>% as.numeric(),
                dmrs$`Clinical Outcome`[index] %>% str_sub(., 2, 1000) %>%
                  #gsub(., "FatPercWB", "WBfatperc") %>%
                  paste0("outcomedelta_", .,  "_", dmrs$`Time Point`[index], "0"),
                "Weight (kg)", padding = 100, nbreaks = 5
)



dmrs[10, ]
make_bl_dmrplot(dmrs$`DMR Chromosome`[10] %>% paste0("chr", .),
                dmrs$`DMR Start`[10] %>% as.numeric(),
                dmrs$`DMR End`[10] %>% as.numeric(),
                "outcomedelta_weight_30", 200
                )

dmrs[25, ]
make_bl_dmrplot(dmrs$`DMR Chromosome`[25] %>% paste0("chr", .),
                dmrs$`DMR Start`[25] %>% as.numeric(),
                dmrs$`DMR End`[25] %>% as.numeric(),
                "outcomedelta_leptin_60", 200
)


dmrs[26, ]
dmrs[26, "Top DMP"]
make_bl_dmrplot(dmrs$`DMR Chromosome`[26] %>% paste0("chr", .),
                dmrs$`DMR Start`[26] %>% as.numeric(),
                dmrs$`DMR End`[26] %>% as.numeric(),
                "outcomedelta_tg_30", 200
)


dmrs[28, ]
make_bl_dmrplot(dmrs$`DMR Chromosome`[28] %>% paste0("chr", .),
                dmrs$`DMR Start`[28] %>% as.numeric(),
                dmrs$`DMR End`[28] %>% as.numeric(),
                "outcomedelta_crprotein_60", 20
)



dmrs[29, ]
make_bl_dmrplot(dmrs$`DMR Chromosome`[29] %>% paste0("chr", .),
                dmrs$`DMR Start`[29] %>% as.numeric(),
                dmrs$`DMR End`[29] %>% as.numeric(),
                "outcomedelta_ldl_60", 200
)





# alt representation as heatmap ------------------------------------------------

library(ComplexHeatmap)

make_bl_dmrplot <- 
  function(dmr_chr, dmr_left, dmr_right, outcome_in, padding,
           custom_outcome_name = "Outcome",
           xlab = paste0(dmr_chr, ":", dmr_left, "-", dmr_right, "\n(length = ", dmr_right - dmr_left, ")" ),
           ylab = "Fraction Methylated (β)", nbreaks = 5) {
    
    probe_in_df <-
      EPIC.hg38.manifest_df %>% 
      filter(seqnames == dmr_chr, between(start, dmr_left - padding, dmr_right  + padding)) %>%
      filter(probe %in% probenames_filtered) %>%
      transmute(probe, pos = start, indmr =  between(pos, dmr_left, dmr_right))
    
    categ_outcome <-
      metadata_final_sample[ , outcome_in, drop = T] %>%
      cut_interval(., n = nbreaks) %>%
      `levels<-`(c(1:nbreaks))
    
    metadata_final_sample$categ_outcome <- categ_outcome
    
    probe_val <-
      zvals_ruvg_cellcomp[ , probenames_filtered %in% probe_in_df$probe]
    # methyl_betas_filtered[ probenames_filtered %in% probe_in_df$probe, ] %>%
    #             as_tsibble %>% .[ , metadata_final_sample$Sample_Name] %>%
    # t
    # betas_raw[ probenamesall %in% probe_in_df$probe, ] %>%
    # as_tibble %>% .[ , metadata_final_sample$Sample_Name] %>%
    # t
    
    tmp_plot <-
      metadata_final_sample %>%
      select(IID, categ_outcome, Timepoint) %>%
      bind_cols(., probe_val %>% as_tibble %>% set_names(probe_in_df$probe)) %>%
      # filter(categ_outcome %in% c(1, 4)) %>%
      # group_by(IID, categ_outcome, Timepoint) %>%
      group_by(categ_outcome, Timepoint) %>%
      #summarize_if(is.numeric, function(x) { median(x, na.rm = T) }) %>%
      pivot_longer(cols = names(.) %>% .[grepl("cg", .)], names_to = "probe") %>%
      left_join(probe_in_df, by = "probe")
    
    heatmap_tmp <- metadata_final_sample %>% filter(Timepoint == "BL") %>%
      .[ , outcome_in] %>% as.data.frame %>% set_names("delta")
    
    sortorder <- heatmap_tmp$delta %>% order
    
    heatmap_tmp <- heatmap_tmp %>%
      mutate(delta = cut_number(delta, n = 4))
    
    Heatmap(matrix = probe_val[metadata_final_sample$Timepoint == "BL", ] %>% .[sortorder, ],
            cluster_rows = T, cluster_columns = F, row_split = heatmap_tmp[sortorder, ]$delta,
            col = circlize::colorRamp2(breaks = c(quantile(probe_val, 0.05), quantile(probe_val, 0.95)), c("white", "black"))) +
      rowAnnotation(df =  heatmap_tmp[sortorder, ])
    
  }

index <- 5
dmrs[index, ]
make_bl_dmrplot(dmrs$`DMR Chromosome`[index] %>% paste0("chr", .),
                dmrs$`DMR Start`[index] %>% as.numeric(),
                dmrs$`DMR End`[index] %>% as.numeric(),
                dmrs$`Clinical Outcome`[index] %>% str_sub(., 2, 1000) %>%
                  #gsub(., "FatPercWB", "WBfatperc") %>%
                  paste0("outcomedelta_", .,  "_", dmrs$`Time Point`[index], "0"),
                "Weight (kg)", padding = 100, nbreaks = 4
)

