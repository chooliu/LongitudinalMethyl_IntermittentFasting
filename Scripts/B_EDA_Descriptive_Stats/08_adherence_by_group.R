# ==============================================================================
# B/08_adherence_by_group.R
# check self-reported adherence levels
# ==============================================================================




# load data, libraries ---------------------------------------------------------

# Q3: dietary adherence level during the past 30 days.
# Q4: how hard it was to adhere to your prescribed study diet during the past 30 days.
# Q5: can adhere to your prescribed diet for the next 30 days.

filepath_adherence <-
  "./Data/20210209-Adherence/DRIFT2Database-BorengasserAdherence_DATA_2021-02-08_1157.csv"

adherence <-
  read_csv(filepath_adherence) %>%
  mutate(time = redcap_event_name %>% as_factor() %>%
           fct_recode(
             `1` = "week_4_arm_1", `2` = "week_8_arm_1", # adherence-only times
             `3` = "week_13_arm_1", `6` = "week_26_arm_1") %>%
           as.character %>% as.numeric,
         record_id = str_pad(record_id, width = 3, side = "left", pad = "0")) %>%
  filter(!is.na(time)) %>%
  pivot_wider(id_cols = record_id,
              names_from = time,
              names_prefix = "month_",
              values_from =
                c("imf_3", "dcr_3", "imf_4", "dcr_4", "imf_5", "dcr_5"))

adherence_final <-
  adherence %>%
  rowwise() %>%
  transmute(record_id, 
            adherence_q3_m1 = mean(c(imf_3_month_1, dcr_3_month_1), na.rm = T),
            adherence_q3_m2 = mean(c(imf_3_month_2, dcr_3_month_2), na.rm = T),
            adherence_q3_m3 = mean(c(imf_3_month_3, dcr_3_month_3), na.rm = T),
            adherence_q3_m6 = mean(c(imf_3_month_6, dcr_3_month_6), na.rm = T),
            adherence_q4_m1 = mean(c(imf_4_month_1, dcr_4_month_1), na.rm = T),
            adherence_q4_m2 = mean(c(imf_4_month_2, dcr_4_month_2), na.rm = T),
            adherence_q4_m3 = mean(c(imf_4_month_3, dcr_4_month_3), na.rm = T),
            adherence_q4_m6 = mean(c(imf_4_month_6, dcr_4_month_6), na.rm = T),
            adherence_q5_m1 = mean(c(imf_5_month_1, dcr_5_month_1), na.rm = T),
            adherence_q5_m2 = mean(c(imf_5_month_2, dcr_5_month_2), na.rm = T),
            adherence_q5_m3 = mean(c(imf_5_month_3, dcr_5_month_3), na.rm = T),
            adherence_q5_m6 = mean(c(imf_5_month_6, dcr_5_month_6), na.rm = T),
            
            meanadherence_q3 = mean(c(imf_3_month_1, imf_3_month_2, imf_3_month_3,
                                      dcr_3_month_1, dcr_3_month_2, dcr_3_month_3),
                                    na.rm = T),
            meanadherence_q4 = mean(c(imf_4_month_1, imf_4_month_2, imf_4_month_3,
                                      dcr_4_month_1, dcr_4_month_2, dcr_4_month_3),
                                    na.rm = T),
            meanadherence_q5 = mean(c(imf_5_month_1, imf_5_month_2, imf_5_month_3,
                                      dcr_5_month_1, dcr_5_month_2, dcr_5_month_3),
                                    na.rm = T)
  ) %>%
  inner_join(metadata_final_subject, by = "record_id")



# adherence at month 1 ---------------------------------------------------------

ggplot(adherence_final,
       aes(as.factor(Intervention), adherence_q3_m1)) +
  geom_boxplot() +
  theme_few()

ggplot(adherence_final,
       aes(as.factor(Intervention), adherence_q4_m1)) +
  geom_boxplot() +
  theme_few()

ggplot(adherence_final,
       aes(as.factor(Intervention), adherence_q5_m1)) +
  geom_boxplot() +
  theme_few()



# adherence at month 3 ---------------------------------------------------------

ggplot(adherence_final,
       aes(as.factor(Intervention), adherence_q3_m2)) +
  geom_boxplot() +
  theme_few()

ggplot(adherence_final,
       aes(as.factor(Intervention), adherence_q4_m2)) +
  geom_boxplot() +
  theme_few()

ggplot(adherence_final,
       aes(as.factor(Intervention), adherence_q5_m2)) +
  geom_boxplot() +
  theme_few()


# over entire trial ------------------------------------------------------------

ggplot(adherence_final,
       aes(as.factor(Intervention), meanadherence_q3)) +
  geom_boxplot() +
  theme_few()

ggplot(adherence_final,
       aes(as.factor(Intervention), meanadherence_q4)) +
  geom_boxplot() +
  theme_few()

ggplot(adherence_final,
       aes(as.factor(Intervention), meanadherence_q5)) +
  geom_boxplot() +
  theme_few()

