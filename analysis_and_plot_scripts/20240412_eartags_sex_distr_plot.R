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
               sessioninfo,
               tidyverse)
# Session info ----
platform_info() %>%
  as.data.frame() %>%
  gather(setting, value, version, os, system, ui, language, collate, ctype, tz, date, rstudio, pandoc) %>%
  add_column("#"="#", .before = "setting") %>%
  print(., row.names = FALSE, right = FALSE)
# setting  value                                                  
# version  R version 4.3.2 (2023-10-31)                           
# os       Ubuntu 22.04.3 LTS                                     
# system   x86_64, linux-gnu                                      
# ui       RStudio                                                
# language en_CA:en                                               
# collate  en_CA.UTF-8                                            
# ctype    en_CA.UTF-8                                            
# tz       America/Vancouver                                      
# date     2023-11-09                                             
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)             
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
  mutate(age = gsub("newborn",
                    0,
                    age_of_calf_at_collection_rounded_weeks)) %>%
  mutate(age = as.numeric(gsub("unknown",
                               NA,
                               age)))
older_calves <- sample_list_single_pathologies %>%
  filter(age > 4) %>%
  select(et_index,
         age) %>%
  mutate(older_calf_outliers = et_index) %>%
  mutate(occurences_of_ages = c(3.5,2,1,1,1,2.5))
sample_list_single_pathologies <- left_join(sample_list_single_pathologies,
                                              older_calves,
                                              by = c("et_index",
                                                     "age"))
unique(sample_list_single_pathologies$pathology)
# [1] "control"    "large_calf" "unknown"    "abnormal" 
sex_dist_of_Abn <- sample_list_single_pathologies %>%
  filter(pathology == "abnormal") %>%
  count(sex)
sex_dist_of_cont <- sample_list_single_pathologies %>%
  filter(pathology == "control") %>%
  count(sex)
sex_dist_of_LOS <- sample_list_single_pathologies %>%
  filter(pathology == "large_calf") %>%
  count(sex)
sex_dist_of_unk <- sample_list_single_pathologies %>%
  filter(pathology == "unknown") %>%
  count(sex)
order_of_pathologies <- c("Control\nCalves",
                          "Large Offspring\nSyndrome",
                          "Morphological\nAbnormalities",
                          "Condition\nNot Documented")
# density plot ----
densplot <- ggplot(data = sample_list_single_pathologies,
                   aes(x = factor(pathology,
                                  levels = c("control",
                                             "large_calf",
                                             "abnormal",
                                             "unknown")),
                       fill = sex)) +
  scale_x_discrete(name = "Calf Conditions",
                   labels = order_of_pathologies) +
  scale_y_continuous(name = "Number of Calves",
                     expand = expansion(),
                     limits = c(0,24),
                     breaks = seq(0,24,by = 2)) +
  geom_bar() +
  scale_fill_manual(values = c("#209e2c",
                               "#841ab5",
                               "grey40"),
                    name = "Calf Sex",
                    labels = c("Female",
                               "Male",
                               "Sex Not Recorded")) +
  geom_text(aes(x = 1,
                y = 13,
                label = sex_dist_of_cont$n[1]), # Control females
            na.rm = TRUE,
            size = 9,
            colour = "#FFFFFF") +
  geom_text(aes(x = 1,
                y = 3,
                label = sex_dist_of_cont$n[2]), # Control males
            na.rm = TRUE,
            size = 9,
            colour = "#FFFFFF") +
  geom_text(aes(x = 2,
                y = 15,
                label = sex_dist_of_LOS$n[1]), #LOS females
            na.rm = TRUE,
            size = 9,
            colour = "#FFFFFF") +
  geom_text(aes(x = 2,
                y = 4.5,
                label = sex_dist_of_LOS$n[2]), #LOS males
            na.rm = TRUE,
            size = 9,
            colour = "#FFFFFF") +
  geom_text(aes(x = 2,
                y = 0.5,
                label = sex_dist_of_LOS$n[3]), #LOS unknowns
            na.rm = TRUE,
            size = 6,
            colour = "#FFFFFF") +
  geom_text(aes(x = 3,
                y = 10.5,
                label = sex_dist_of_Abn$n[1]), #Abnormal females
            na.rm = TRUE,
            size = 9,
            colour = "#FFFFFF") +
  geom_text(aes(x = 3,
                y = 4,
                label = sex_dist_of_Abn$n[2]), #Abnormal males
            na.rm = TRUE,
            size = 9,
            colour = "#FFFFFF") +
  geom_text(aes(x = 3,
                y = 0.5,
                label = sex_dist_of_Abn$n[3]), #Abnormal unknowns
            na.rm = TRUE,
            size = 6,
            colour = "#FFFFFF") +
  geom_text(aes(x = 4,
                y = 8,
                label = sex_dist_of_unk$n[1]), #unknown females
            na.rm = TRUE,
            size = 9,
            colour = "#FFFFFF") +
  geom_text(aes(x = 4,
                y = 4.5,
                label = sex_dist_of_unk$n[2]), #unknown males
            na.rm = TRUE,
            size = 9,
            colour = "#FFFFFF") +
  geom_text(aes(x = 4,
                y = 1,
                label = sex_dist_of_unk$n[3]), #unknown unknowns
            na.rm = TRUE,
            size = 9,
            colour = "#FFFFFF") +
  theme_classic() +
  # ggtitle("Distribution of Calf Sex by Condition") +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16,
                                  colour = "grey20"),
        axis.title.x = element_text(face = "bold",
                                    size = 14,
                                    colour = "grey20"),
        axis.text.x = element_text(face = "bold",
                                   size = 16,
                                   colour = "grey30"),
        axis.title.y = element_text(face = "bold",
                                    size = 14,
                                    colour = "grey20"),
        axis.text.y = element_text(face = "bold",
                                   size = 16,
                                   colour = "grey30"),
        legend.title = element_text(face = "bold",
                                   size = 16,
                                   colour = "grey20"),
        legend.text = element_text(face = "bold",
                                   size = 18,
                                   colour = "grey30"),
        legend.position = c(0.87,0.83),
        plot.margin = margin(0.5,1,0.5,0.5, "cm")
        )
densplot  
ggsave("data/analyzed_data/20240412_distrb_of_calves_sex_grouped_conditions.png",
       dpi = 300,
       plot = densplot)
