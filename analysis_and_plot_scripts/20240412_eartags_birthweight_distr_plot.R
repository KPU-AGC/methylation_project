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
  drop_na(calf_birth_weight_lbs)
# Number of calves with unknown weights ----
number_of_unknown_weights <- as.data.frame(table(sample_list_single_pathologies$calf_birth_weight_lbs == "unknown")) %>%
  filter(Var1 == "TRUE") %>%
  select(Freq) %>%
  as.numeric()
# Data frame of only calves with known weights ----
calves_w_known_birth_weight <- sample_list_single_pathologies %>%
  filter(!calf_birth_weight_lbs == "unknown") %>%
  mutate_at("calf_birth_weight_lbs",
            as.numeric) # Had to transform the column to numeric as the csv is imported with each column being character in this script for some reason
str(calves_w_known_birth_weight)
# Data frame with just the calves that have unknown weights ---
calves_w_unknown_birth_weight <- sample_list_single_pathologies %>%
  filter(calf_birth_weight_lbs == "unknown")
# What's the distribution of calves' conditions with unknown weights? ----
dist_of_calves_w_unknown_weights <- as.data.frame(table(calves_w_unknown_birth_weight$pathology))
# Exactly what are the conditions/pathologies labelled as in the data frame? ----
unique(sample_list_single_pathologies$pathology)
# [1] "control"    "large_calf" "unknown"    "abnormal"
# density plot ----
densplot <- ggplot(data = calves_w_known_birth_weight,
                   aes(x = calf_birth_weight_lbs,
                       fill = pathology)) +
  scale_x_continuous(name = "Birth Weight of Calves (lbs)",
                     expand = expansion(),
                     limits = c(10,170),
                     breaks = seq(10,170,by = 10)) +
  scale_y_continuous(name = "Number of Calves",
                     expand = expansion(),
                     limits = c(0,8),
                     breaks = seq(0,8,by = 2)) +
  geom_histogram(binwidth = 3) +
  # geom_vline(xintercept = 90,
  #            linetype = 3,
  #            linewidth = 1.5,
  #            colour = "#e65603") +
  # geom_text(aes(x = 87,
  #               y = 4,
  #               label = "Average Weight"),
  #           na.rm = TRUE,
  #           size = 6,
  #           colour = "#e65603",
  #           angle = 90) +
  scale_fill_manual(values = c("#840029",
                                 "#006699",
                                 "#f2a40d",
                                 "grey40"),
                    name = "Condition",
                    labels = c("Abnormal Calves",
                               "Control Calves",
                               "Large Offspring (LOS)",
                               "Not Documented")) +
  theme_classic() +
  # ggtitle("Distribution of Calf Ages at Time of Sample Collection") +
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
        legend.position = c(0.75,0.75),
        plot.margin = margin(0.25,0.5,0,0,"cm")
        # plot.margin = margin(0.5,1,0.5,0.5, "cm")
        )
densplot  
ggsave("data/analyzed_data/20240412_distrb_of_calves_weight_grouped_conditions_sans_avg_weight.png",
       dpi = 300,
       plot = densplot)
