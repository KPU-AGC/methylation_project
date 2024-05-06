# ----
#                    ▓▓░  ▓▓▓  ░▓▓▓
#   -. .-.   .-. .- ▓░ ▓░▓░   ░▓  ░▓ -. .-.   .-. .-
#   ||\|||\ /|||\|| ▓▓▓▓░▓░░▓▓░▓     ||\|||\ /|||\||
#   |/ \|||\|||/ \| ▓░ ▓░▓░ ░▓░▓  ░▓ |/ \|||\|||/ \|
#   ¯   `-¯ `-'   ` ▓░ ▓░ ▓▓▓  ░▓▓▓  ¯   `-¯ `-'   `
#  T H E  A P P L I E D  G E N O M I C S  C E N T R E
#           --- Coding a better harvest ---
# Background ----
# This script is for plotting out the distributions of ages of calves that were sampled for our ear tag NGS data
# The purpose of this script is to have an idea of the range of ages
# And therefore provide support for the rationale that it will be acceptable to parse data based on inclusion of just "newborns"/neonatal calves
# A definite age for which cattle are considered neonates doesn't appear to be out there from the perspective of regulations and handling of livestock
# However, an age of 1 week or less aligns with Canadian regulations in which animals are universally still considered neonates
# For this analysis here I'm going to stick with using calves that are less than a week old as neonates
# Loading packages ----
library(pacman)
pacman::p_load(imgpalr,
               methylKit,
               patchwork,
               reshape2,
               sessioninfo,
               tidyverse)
# Session info ----
platform_info() %>%
  as.data.frame() %>%
  gather(setting, value, version, os, system, ui, language, collate, ctype, tz, date, rstudio, pandoc) %>%
  add_column("#"="#", .before = "setting") %>%
  print(., row.names = FALSE, right = FALSE)
# setting  value                                                  
# version  R version 4.3.3 (2024-02-29)                           
# os       Ubuntu 22.04.4 LTS                                     
# system   x86_64, linux-gnu                                      
# ui       RStudio                                                
# language en_CA:en                                               
# collate  en_CA.UTF-8                                            
# ctype    en_CA.UTF-8                                            
# tz       America/Vancouver                                      
# date     2024-04-11                                             
# rstudio  2023.09.1+494 Desert Sunflower (desktop)               
# pandoc   3.1.3 @ /home/patrick/miniconda3/envs/R4.3.0/bin/pandoc
set.seed(123)
rnorm(5) %>%
  as.data.frame() %>%
  rename(., "rnorm_values" = ".") %>%
  add_column("#"="#",
             .before = "rnorm_values") %>%
  print(.,
        row.names = FALSE,
        right = FALSE)
# rnorm_values
# -0.56047565 
# -0.23017749 
#  1.55870831 
#  0.07050839 
#  0.12928774
# Loading metadata table which contains the age ranges ----
sample_list_single_pathologies <- read_csv("data/analyzed_data/metadata_formatted_for_parsing_single_pathologies.csv") %>%
  mutate(age = gsub("newborn", # Newborn or neonate was taken as less than 1 week, i.e., 0 weeks old, based on what I could find from the Government of Canada's legislation on the topic: https://laws-lois.justice.gc.ca/eng/regulations/C.R.C.,_c._296/index.html I have more detailed description of this in my notes from 2023-11-09 in my 2023 labbook
                    0,
                    age_of_calf_at_collection_rounded_weeks)) %>%
  mutate(age = as.numeric(gsub("unknown",
                               NA,
                               age))) %>%
  drop_na(age)
older_calves <- sample_list_single_pathologies %>%
  filter(age > 4) %>% # Age of 4 weeks was taken simply because less than this age is where the majority of samples lie 
  select(et_index,
         age) %>%
  mutate(older_calf_outliers = et_index) %>%
  mutate(occurences_of_ages = c(3.5,2,1,1,1,2.5))
sample_list_single_pathologies <- left_join(sample_list_single_pathologies,
                                              older_calves,
                                              by = c("et_index",
                                                     "age"))
# Creation of dataframe for stacked histogram ----
cont_dist <- sample_list_single_pathologies %>%
  filter(pathology == "control")
table(cont_dist$age)
LOS_dist <- sample_list_single_pathologies %>%
  filter(pathology == "large_calf")
table(LOS_dist$age)
abn_dist <- sample_list_single_pathologies %>%
  filter(pathology == "abnormal")
table(abn_dist$age)
unk_dist <- sample_list_single_pathologies %>%
  filter(pathology == "unknown")
table(unk_dist$age)

age_dist_df <- data.frame(age_range = c(min(sample_list_single_pathologies$age):max(sample_list_single_pathologies$age)),
                          cont_dist_of_ages = c(14, 2, 1, 0, 0, 0, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                1, 0),
                          LOS_dist_of_ages = c(18, 2, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0),
                          abn_dist_of_ages = c(10, 0, 2, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 1, 0, 0, 0, 0, 0, 1,
                                               0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 1),
                          unknown_of_ages = c(5, 1, 0, 0, 1, 0, 0, 0, 0, 0,
                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                              0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                              0, 0))
# Stacked histogram plot ----
stacked_histogram <- ggplot(data = age_dist_df,
                            aes(x = age_range,
                                fill = L1)) +
    geom_histogram(position = "stack") #+
stacked_histogram  
ggsave("data/analyzed_data/distrb_of_calf_ages_with_older_calves_boxed.png",
       plot = densplot)
