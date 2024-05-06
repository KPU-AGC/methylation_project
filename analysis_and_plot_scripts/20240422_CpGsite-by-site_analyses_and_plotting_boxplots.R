# ----
#                    ▓▓░  ▓▓▓  ░▓▓▓
#   -. .-.   .-. .- ▓░ ▓░▓░   ░▓  ░▓ -. .-.   .-. .-
#   ||\|||\ /|||\|| ▓▓▓▓░▓░░▓▓░▓     ||\|||\ /|||\||
#   |/ \|||\|||/ \| ▓░ ▓░▓░ ░▓░▓  ░▓ |/ \|||\|||/ \|
#   ¯   `-¯ `-'   ` ▓░ ▓░ ▓▓▓  ░▓▓▓  ¯   `-¯ `-'   `
#  T H E  A P P L I E D  G E N O M I C S  C E N T R E
#           --- Coding a better harvest ---
# Background ----
# While not critical to our analyses as the violin plots of mean percent methylation for each sample seems to give the right level of focus on our data
# This here is a script for performing PCA analyses where we can look at the global methylation for each sample
# The purpose of this script is to give us the resolution of whole conditions asking the experimental question if there are clusterings of conditions 
# Loading packages ----
library(pacman)
pacman::p_load(ggh4x,
               patchwork,
               rstatix,
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
# date     2024-04-24                                             
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
# Loading functions ----
# Source functions ----
source("~/Documents/general_code/getOutliers_function.R")
source("code/function_generate-mean-total-count-CpG-sites.R")
source("code/function_pull-all-mean-total-count-together.R")
source("code/function_filter-low-rep-CpG-sites.R")
source("code/function_pull-all-low-rep-CpG-sites-together.R")
# Loading Metadata ----
sample_list <- read_csv("data/analyzed_data/metadata_formatted_for_parsing_single_pathologies.csv")
target_list <- read_csv("data/analyzed_data/bsx_target_sites_trimmed_target_name.csv") %>%
  rename("target_site" = "primer")
# Loading data ----
full_run4_dataset <- read_tsv("data/raw_data/20230920_195755_nonCG_filtered/all_samples_together_S1_L001_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.mod_xsomal_loci.all_insert_target_name.cov",
                              col_names = FALSE) %>%
  rename("et_index" = "X1",
         "target_site" = "X2",
         "CpG_site" = "X3",
         "perc_meth" = "X4",
         "meth_C_count" = "X5",
         "unmeth_C_count" = "X6") %>%
  mutate(total_count = meth_C_count + unmeth_C_count) %>%
filter(!grepl("NC_",
              target_site)) %>%
filter(!grepl("NW_",
              target_site))
off_target_sites_run4_dataset <- read_tsv("data/raw_data/20230920_195755_nonCG_filtered/all_samples_together_S1_L001_1_bismark_bt2_pe.nonCG_filtered.bedGraph.gz.bismark.zero.mod_xsomal_loci.all_insert_target_name.cov",
                                          col_names = FALSE) %>%
  rename("et_index" = "X1",
         "target_site" = "X2",
         "CpG_site" = "X3",
         "perc_meth" = "X4",
         "meth_C_count" = "X5",
         "unmeth_C_count" = "X6") %>%
  filter(grepl("_", target_site)) %>%
  mutate(total_count = meth_C_count + unmeth_C_count) %>%
  left_join(.,
            sample_list,
            by = "et_index")
# Calculate means of total_counts for each target site ----
pull_all_mean_total_count_together(full_run4_dataset)
pull_all_mean_total_count_together(off_target_sites_run4_dataset) # Adding in this one to calculate some metrics of the off-target sites
# Calculate some metrics for off-target sites ----
# Things like number of sites represented,
n_distinct(all_mean_total_count_together_offTarg$target_site)
# [1] 5877
# average read depth,
mean(all_mean_total_count_together_offTarg$total_count)
# [1] 12.57299
# and total portion of sequencing reads from the run
sum(all_mean_total_count_together_offTarg$total_count)
# [1] 208171
# Compare this to the on-target sites ----
n_distinct(all_mean_total_count_together$target_site)
# [1] 20
# Good. This is what we should expect
mean(all_mean_total_count_together$total_count)
# [1] 9779.869
sum(all_mean_total_count_together$total_count)
# [1] 153064737
paste0(((sum(all_mean_total_count_together_offTarg$total_count))/(sum(full_run4_dataset$total_count)))*100,"%")
# [1] "0.134304412517792%"
# The off-target sites represent just 0.13% of CpG reads from the whole sequencing run
# Parse out total_counts that are less than 25% the mean value/representation of CpG sites in the target site ----
pull_all_low_rep_CpG_sites_analyses_together(all_mean_total_count_together)
only_hi_rep_CpG_sites <- all_low_rep_CpG_sites_analyses_together %>%
  filter(low_rep_CpG_site == FALSE) %>%
  filter(total_count > 15) 
  # select(-c(mean_total_count,
  #           min_total_count,
  #           max_total_count,
  #           range_total_count))
only_lo_rep_CpG_sites <- all_low_rep_CpG_sites_analyses_together %>%
  filter(low_rep_CpG_site == TRUE) %>%
  filter(total_count > 15) #%>%
  # select(-c(mean_total_count,
  #           min_total_count,
  #           max_total_count,
  #           range_total_count))
only_hi_rep_recalc_tot_count_CpG_sites <- pull_all_mean_total_count_together(only_hi_rep_CpG_sites)
only_lo_rep_recalc_tot_count_CpG_sites <- pull_all_mean_total_count_together(only_lo_rep_CpG_sites)
# write_csv("data/analyzed_data/20240423_CpG_sites_results_w_low_rep_removed.csv")
only_hi_rep_recalc_tot_count_CpG_sites_for_plots <- left_join(only_hi_rep_recalc_tot_count_CpG_sites,
                                                    sample_list,
                                                    by = "et_index") %>%
  filter(!pathology %in% c("abnormal","unknown"))
only_hi_rep_recalc_tot_count_CpG_sites_for_plots$CpG_site <- as.character(only_hi_rep_recalc_tot_count_CpG_sites_for_plots$CpG_site)
class(only_hi_rep_recalc_tot_count_CpG_sites_for_plots$CpG_site)
only_lo_rep_recalc_tot_count_CpG_sites_for_plots <- left_join(only_lo_rep_recalc_tot_count_CpG_sites,
                                                              sample_list,
                                                              by = "et_index") %>%
  filter(!pathology %in% c("abnormal","unknown"))
only_lo_rep_recalc_tot_count_CpG_sites_for_plots$CpG_site <- as.character(only_lo_rep_recalc_tot_count_CpG_sites_for_plots$CpG_site)
class(only_lo_rep_recalc_tot_count_CpG_sites_for_plots$CpG_site)
all_low_rep_CpG_sites_analyses_together_for_plots <- left_join(all_low_rep_CpG_sites_analyses_together,
                                                     sample_list,
                                                     by = "et_index") %>%
  filter(!pathology %in% c("abnormal","unknown"))
all_low_rep_CpG_sites_analyses_together_for_plots$CpG_site <- as.character(all_low_rep_CpG_sites_analyses_together_for_plots$CpG_site)
class(all_low_rep_CpG_sites_analyses_together_for_plots$CpG_site)
# Check to see if there's a disparity between filtered and low representation CpG sites ----
# (At a glance it looks like there's more LOS calves in the low represented sites)
number_of_CpG_sites_ctrl_filt_data <- only_hi_rep_recalc_tot_count_CpG_sites_for_plots %>%
  filter(pathology == "control")
paste0("There are ",length(unique(number_of_CpG_sites_ctrl_filt_data$CpG_site))," unique CpG sites represented among control/healthy calves in filtered data")
# [1] "There are 332 unique CpG sites represented among control/healthy calves in filtered data"
number_of_CpG_sites_LOS_filt_data <- only_hi_rep_recalc_tot_count_CpG_sites_for_plots %>%
  filter(pathology == "large_calf")
paste0("There are ",length(unique(number_of_CpG_sites_LOS_filt_data$CpG_site))," unique CpG sites represented among LOS calves in filtered data")
# [1] "There are 333 unique CpG sites represented among LOS calves in filtered data"
number_of_CpG_sites_ctrl_lo_rep_data <- only_lo_rep_recalc_tot_count_CpG_sites_for_plots %>%
  filter(pathology == "control")
paste0("There are ",length(unique(number_of_CpG_sites_ctrl_lo_rep_data$CpG_site))," unique CpG sites represented among control/healthy calves in low represented CpG site data")
# [1] "There are 335 unique CpG sites represented among control/healthy calves in low represented CpG site data"
number_of_CpG_sites_LOS_lo_rep_data <- only_lo_rep_recalc_tot_count_CpG_sites_for_plots %>%
  filter(pathology == "large_calf")
paste0("There are ",length(unique(number_of_CpG_sites_LOS_lo_rep_data$CpG_site))," unique CpG sites represented among LOS calves in low represented CpG site data")
# [1] "There are 333 unique CpG sites represented among LOS calves in low represented CpG site data"
# What's this extra CpG site in the LOS calves? Is it the same in both these data sets?
unique_CpG_sites_ctrl_filt_data <- data.frame(unique(number_of_CpG_sites_ctrl_filt_data$CpG_site)) %>%
  rename(.,
         CpG_site = unique.number_of_CpG_sites_ctrl_filt_data.CpG_site.)
unique_CpG_sites_LOS_filt_data <- data.frame(unique(number_of_CpG_sites_LOS_filt_data$CpG_site)) %>%
  rename(.,
         CpG_site = unique.number_of_CpG_sites_LOS_filt_data.CpG_site.)
LOS_unique_site_filt_data <- unique_CpG_sites_LOS_filt_data %>%
  anti_join(.,
            unique_CpG_sites_ctrl_filt_data,
            by = "CpG_site")
ctrl_unique_site_filt_data <- unique_CpG_sites_ctrl_filt_data %>%
  anti_join(.,
            unique_CpG_sites_LOS_filt_data,
            by = "CpG_site")
unique_CpG_sites_ctrl_lo_rep_data <- data.frame(unique(number_of_CpG_sites_ctrl_lo_rep_data$CpG_site)) %>%
  rename(.,
         CpG_site = unique.number_of_CpG_sites_ctrl_lo_rep_data.CpG_site.)
unique_CpG_sites_LOS_lo_rep_data <- data.frame(unique(number_of_CpG_sites_LOS_lo_rep_data$CpG_site)) %>%
  rename(.,
         CpG_site = unique.number_of_CpG_sites_LOS_lo_rep_data.CpG_site.)
LOS_unique_site_lo_rep <- unique_CpG_sites_LOS_lo_rep_data %>%
  anti_join(.,
            unique_CpG_sites_ctrl_lo_rep_data,
            by = "CpG_site")
ctrl_unique_site_lo_rep <- unique_CpG_sites_ctrl_lo_rep_data %>%
  anti_join(.,
            unique_CpG_sites_LOS_lo_rep_data,
            by = "CpG_site")
paste0("There's 1 CpG site found in IGF2R in the LOS calves that aren't in the controls for the filtered data: ",list(LOS_unique_site_filt_data$CpG_site))
# [1] "There's 1 CpG site found in IGF2R in the LOS calves that aren't in the controls for the filtered data: 96223335"
paste0("It boils down to the fact that there's 4 CpG sites found in LIF in the controls that aren't in the LOS calves for the low represented data: ",list(ctrl_unique_site_lo_rep$CpG_site))
paste0("And 2 sites in RTL1 in the LOS calves not in the controls: for the low represented data: ",list(LOS_unique_site_lo_rep$CpG_site))
# [1] "It boils down to the fact that there's 4 CpG sites found in LIF  in the controls that aren't in the LOS calves for the low represented data: 69264124, 69264133, 69264295, 65778664
# [1] "And 2 sites in IGF2R and MEST, respectively, in the LOS calves not in the controls: for the low represented data: 96223294, 94250168"
paste0("There are ",(length(which(only_hi_rep_recalc_tot_count_CpG_sites_for_plots$pathology == "control")))," CpG sites total from control/healthy calves in the filtered data of CpG sites")
# [1] "There are 4548 CpG sites total from control/healthy calves in the filtered data of CpG sites"
paste0("There are ",(length(which(only_hi_rep_recalc_tot_count_CpG_sites_for_plots$pathology == "large_calf")))," CpG sites total from LOS calves in the filtered data of CpG sites")
# [1] "There are 4644 CpG sites total from LOS calves in the filtered data of CpG sites"
paste0("This matches the total number of ",((length(which(only_hi_rep_recalc_tot_count_CpG_sites_for_plots$pathology == "control")))+(length(which(only_hi_rep_recalc_tot_count_CpG_sites_for_plots$pathology == "large_calf"))))," CpG sites represented in this dataset")
# [1] "This matches the total number of 9192 CpG sites represented in this dataset"
paste0("There are ",(length(which(only_lo_rep_recalc_tot_count_CpG_sites_for_plots$pathology == "control")))," CpG site counts from control/healthy calves in the low represented CpG site data")
# [1] "There are 1186 CpG site counts from control/healthy calves in the low represented CpG site data"
paste0("There are ",(length(which(only_lo_rep_recalc_tot_count_CpG_sites_for_plots$pathology == "large_calf")))," CpG site counts from LOS calves in the low represented CpG site data")
# [1] "There are 1556 CpG site counts from LOS calves in the low represented CpG site data"
paste0("This matches the total number of ",((length(which(only_lo_rep_recalc_tot_count_CpG_sites_for_plots$pathology == "control")))+(length(which(only_lo_rep_recalc_tot_count_CpG_sites_for_plots$pathology == "large_calf"))))," CpG sites represented in this dataset")
# [1] "This matches the total number of 2742 CpG sites represented in this dataset"
# Universal odds and ends for plots ----
# order_of_sexes <- c("female" = "Female",
#                     "male" = "Male",
#                     "unknown" = "Not Recorded")
# sex_colours <- c("female" = "#209e2c",
#                  "male" = "#841ab5",
#                  "unknown" = "grey20")
order_of_pathologies <- c("control" = "Control",
                          "large_calf" = "Large Calf")
pathology_colours <- c("control" = "#006699",
                       "large_calf" = "#f2a40d")
# Violin plots for percent methylation of each CpG site ----
for (i in unique(only_lo_rep_recalc_tot_count_CpG_sites_for_plots$target_site)) {
  # Select the specific target and determine outliers ----
  select_target_to_plot <- only_lo_rep_recalc_tot_count_CpG_sites_for_plots %>%
    filter(target_site == i) %>%
    # mutate(outlier = et_index) %>%
    # mutate(outlier = ifelse(is_outlier(perc_meth),
    #                         et_index,
    #                         as.numeric(NA))) %>%
    # mutate(is_Kune = ifelse(et_index == "LLN-Kune",
    #                         et_index,
    #                         as.numeric(NA))) %>%
    mutate(perc_meth_wo_KUNE = ifelse(et_index == "LLN-Kune",
                                      as.numeric(NA),
                                      perc_meth))
  # select_target_to_plot$outlier[which(is.na(select_target_to_plot$perc_meth))] <- as.numeric(NA)
  # Other iterations of these plots (those for mean % methylation of samples based on conditions) involved a step here to calculate stats
  # I haven't included that here because I feel that this is too cluttered and already filtering steps have been taken to remove erroneous sites
  # Stats here won't be particularly telling of any biological differences and is better used just to observe and get an idea of what we're seeing in the data at this point
  # Label the imprinting statuses ----
  imprinted_status <- target_list %>%
    filter(target_site == i) %>%
    mutate_if(is.character, ~replace_na(.,""))
  target_chromosome_number <- target_list %>%
    filter(target_site == i) %>%
    mutate_if(is.numeric, ~replace_na(.,""))
  # Plot the percent methylation ----
  boxplot <- ggplot(data = select_target_to_plot,
                    aes(x = CpG_site,
                        y = perc_meth_wo_KUNE,
                        fill = factor(pathology,
                                      levels = c("control",
                                                 "large_calf")))) +
    scale_y_continuous(name = "Percent Methylation (%)",
                       limits = c(0,100),
                       breaks = c(0,25, 50, 75, 100)) +
    # sec.axis = sec_axis(~ (./yaxis_range_slope)-yaxis_range_intercept,
    #                     name = "Average Range of Read Depth",
    #                     breaks = seq(0,6,by=1),
    #                     labels = c("1","10","100","1,000","10,000","100,000","1,000,000"))) +
    scale_x_discrete() +
    # scale_x_discrete(limits = c("control",
    #                             "large_calf",
    #                             "abnormal",
    #                             "unknown"),
    #                  breaks = c("control",
    #                             "large_calf",
    #                             "abnormal",
    #                             "unknown"),
    #                  labels = order_of_pathologies) +
    # geom_linerange(aes(x = pathology,
    #                    ymin = yaxis_range_slope*(log10(mean_min_read_count))+yaxis_range_intercept,
  #                    ymax = yaxis_range_slope*(log10(mean_max_read_count))+yaxis_range_intercept),
  #                linewidth = 12,
  #                alpha = 0.6,
  #                colour = "#16a864") +
    geom_boxplot(colour = "#2f3b42",
                 width = 0.5,
                 # lwd = 0.2,
                 show.legend = TRUE,
                 position = position_dodge(0.5),
                 outlier.shape = NA) +
    stat_boxplot(geom = 'errorbar',
                 colour = "#2f3b42",
                 width = 0.25,
                 # lwd = 0.2,
                 show.legend = TRUE,
                 position = position_dodge(0.5)) +
    # stat_summary(fun.data = "mean_cl_boot",
    #              geom = "pointrange",
    #              colour = "#2f3b42",
    #              size = 0.5,
    #              show.legend = TRUE) +
    # geom_jitter(aes(x = factor(pathology,
    #                            levels = c("control",
    #                                       "large_calf",
    #                                       "abnormal",
    #                                       "unknown")),
    #                 y = perc_meth,
    #                 colour = sex),
    #             position = position_jitter(0.1),
    #             size = 1,
    #             show.legend = TRUE) +
    # scale_colour_manual(breaks = names(order_of_sexes),
    #                     values = sex_colours,
    #                     name = "Calf Sex",
    #                     labels = order_of_sexes) +
    scale_fill_manual(breaks = names(order_of_pathologies),
                      values = pathology_colours,
                      name = "Pathology",
                      labels = order_of_pathologies) +
    # facet_wrap(~CpG_site,
    #            strip.position = "bottom",
    #            nrow = 1,
    #            scales = "free_x") +
    # geom_segment(aes(x = 1,
    #                  y = max(perc_meth) + 6,
    #                  xend = 4,
    #                  yend = max(perc_meth) + 6),
    #              colour = "grey20") +
    # geom_segment(aes(x = 1,
    #                  y = max(perc_meth) + 2,
    #                  xend = 1,
    #                  yend = max(perc_meth) + 6),
    #              colour = "grey20") +
    # geom_segment(aes(x = 4,
  #                  y = max(perc_meth) + 2,
  #                  xend = 4,
  #                  yend = max(perc_meth) + 6),
  #              colour = "grey20") +
  # geom_text(aes(x = 2.5,
  #               y = max(perc_meth) + 9,
  #               label = paste0("p = ",test_data$p_value)),
  #           size = 4,
  #           colour = "grey30") +
  theme_classic() +
    xlab(paste0("CpG Location on Chromosome ",target_chromosome_number$chromosome_no)) +
    ggtitle(paste0(i,imprinted_status$imprinted_asterisk)) +
    # geom_text(aes(x = CpG_site,
    #               y = perc_meth,
    #               fill = factor(pathology,
    #                             levels = c("control",
    #                                        "large_calf")),
    #               label = outlier),
    #           na.rm = TRUE,
    #           size = 3,
    #           position = position_nudge(x = 0.1),
    #           colour = "#ea0049",
    #           inherit.aes = TRUE) +
    # geom_text(aes(x = factor(pathology,
    #                          levels = c("control",
    #                                     "large_calf",
    #                                     "abnormal",
    #                                     "unknown")),
    #               y = perc_meth,
    #               label = is_Kune),
    #           na.rm = TRUE,
    #           hjust = -0.1,
    #           size = 3,
    #           position = position_nudge(x = 0.1),
    #           colour = "#ea0049") +
    theme(plot.title = element_text(hjust = 0.5,
                                    face = "bold",
                                    size = 18,
                                    colour = "grey20"),
          axis.title.x = element_text(face = "bold",
                                      size = 14,
                                      colour = "grey20"),
          # axis.text.x = element_blank(),
          axis.text.x = element_text(face = "bold",
                                     size = 12,
                                     angle = 45,
                                     hjust = 0.97,
                                     vjust = 0.95,
                                     colour = "grey30"),
          # axis.ticks.x = element_blank(),
          # axis.line.x = element_blank(),
          axis.title.y = element_text(face = "bold",
                                      size = 14,
                                      colour = "grey20"),
          axis.text.y = element_text(face = "bold",
                                     size = 14,
                                     colour = "grey30"),
          # axis.title.y.right = element_text(face = "bold",
          #                                   size = 12,
          #                                   vjust = 1,
          #                                   colour = "grey20"),
          # axis.text.y.right = element_text(face = "bold",
          #                                  size = 8,
          #                                  colour = "grey30"),
          legend.title = element_text(face = "bold",
                                      size = 14,
                                      colour = "grey20"),
          legend.text = element_text(face = "bold",
                                     size = 14,
                                     colour = "grey30"),
          plot.background = element_rect(fill = "transparent",
                                         colour = NA),
          panel.background = element_rect(fill = "transparent",
                                          colour = NA),
          legend.background = element_rect(fill = "transparent",
                                           colour = NA),
          legend.box.background = element_rect(fill = "transparent",
                                               colour = NA))
          # panel.border = element_blank(),
          # strip.text = element_text(face = "bold",
          #                           size = 14,
          #                           angle = 90,
          #                           colour = "grey20"),
          # strip.background = element_rect(fill = "#b5c3ca",
          #                                 colour = NA))
  export_plot <- boxplot +
    plot_annotation(caption = paste0("Notes:\n• ",imprinted_status$imprinted_asterisk,"Target is ",imprinted_status$which_parent_imprinted," imprinted\n• Outliers from the 0.25 to 0.75 quantiles are labelled")) &
    theme(plot.caption = element_text(hjust = 0.00,
                                      size = 14))
  # Save the plot ----
  # ggsave(paste0("data/analyzed_data/boxplots_20240422_ear_tag_all_samples/20240422_target_",i,"_boxplot_all_CpG_sites_incl.png"),
  #        width = 20,
  #        height = 7.5,
  #        dpi = 300,
  #        plot = export_plot)
  # ggsave(paste0("data/analyzed_data/boxplots_20240422_ear_tag_all_samples/20240422_target_",i,"_boxplot_filt_CpG_sites_only.png"),
  #        width = 20,
  #        height = 7.5,
  #        dpi = 300,
  #        plot = export_plot)
  ggsave(paste0("data/analyzed_data/boxplots_20240422_ear_tag_all_samples/20240422_target_",i,"_boxplot_lo_rep_CpG_sites_only.png"),
         width = 20,
         height = 7.5,
         dpi = 300,
         plot = export_plot)
  # assign(paste0("boxplot_all_CpG_sites_incl_",i),
  #        boxplot)
  # assign(paste0("boxplot_lo_rep_CpG_sites_only_",i),
  #        boxplot)
  # assign(paste0("boxplot_filt_CpG_sites_only_",i),
  #        boxplot)
}
# Combined violin plots into one figure----
boxplot_lo_rep_CpG_sites_only_SNRPN
# layout_design <- "
# AABBCCDDE
# FFGGHHIIZ
# JJKKLLMMZ
# NNOOPPQQZ
# RRSSTTUUZ
# "
layout_design <- "
AABBCC
DDEEFF
GGHHII
JJKKLL
MMNNOO
PPQQRR
SSTTU#
"
# combined_boxplots <- boxplot_all_CpG_sites_incl_DMAP1 + boxplot_all_CpG_sites_incl_DNMT3A + boxplot_all_CpG_sites_incl_DNMT3B + boxplot_all_CpG_sites_incl_GNAS + guide_area() +
#   boxplot_all_CpG_sites_incl_H19 + boxplot_all_CpG_sites_incl_IGF2R + boxplot_all_CpG_sites_incl_KCNQ1 + boxplot_all_CpG_sites_incl_LIFR +
#   boxplot_all_CpG_sites_incl_LIF + boxplot_all_CpG_sites_incl_MEST + boxplot_all_CpG_sites_incl_NNAT + boxplot_all_CpG_sites_incl_PEG3 +
#   boxplot_all_CpG_sites_incl_PEG10 + boxplot_all_CpG_sites_incl_PLAGL1 + boxplot_all_CpG_sites_incl_RTL1 + boxplot_all_CpG_sites_incl_SLC2A8 +
#   boxplot_all_CpG_sites_incl_SNRPN + boxplot_all_CpG_sites_incl_SUV39H1 + boxplot_all_CpG_sites_incl_TXNIP + boxplot_all_CpG_sites_incl_XIST +
#   plot_layout(design = layout_design,
#               guides = "collect"
#               # axis_titles = "collect_y"
#               ) &
#   theme(plot.tag = element_text(face = "bold",
#                                 size = 14,
#                                 colour = "grey20",
#                                 hjust = 0.00),
#         plot.background = element_rect(fill = "transparent",
#                                        colour = NA))
# combined_boxplots <- boxplot_filt_CpG_sites_only_DMAP1 + boxplot_filt_CpG_sites_only_TXNIP + boxplot_filt_CpG_sites_only_MEST +
#   boxplot_filt_CpG_sites_only_PEG10 + boxplot_filt_CpG_sites_only_IGF2R + boxplot_filt_CpG_sites_only_PLAGL1 +
#   boxplot_filt_CpG_sites_only_DNMT3A + boxplot_filt_CpG_sites_only_SLC2A8 + boxplot_filt_CpG_sites_only_DNMT3B +
#   boxplot_filt_CpG_sites_only_GNAS + boxplot_filt_CpG_sites_only_NNAT + boxplot_filt_CpG_sites_only_LIF +
#   boxplot_filt_CpG_sites_only_PEG3 + boxplot_filt_CpG_sites_only_LIFR + boxplot_filt_CpG_sites_only_RTL1 +
#   boxplot_filt_CpG_sites_only_SNRPN + boxplot_filt_CpG_sites_only_H19 + boxplot_filt_CpG_sites_only_KCNQ1 +
#   boxplot_filt_CpG_sites_only_SUV39H1 + boxplot_filt_CpG_sites_only_XIST + guide_area() +
combined_boxplots <- boxplot_lo_rep_CpG_sites_only_DMAP1 + boxplot_lo_rep_CpG_sites_only_DNMT3A + boxplot_lo_rep_CpG_sites_only_DNMT3B +
  boxplot_lo_rep_CpG_sites_only_GNAS + boxplot_lo_rep_CpG_sites_only_H19 + boxplot_lo_rep_CpG_sites_only_IGF2R +
  boxplot_lo_rep_CpG_sites_only_KCNQ1 + boxplot_lo_rep_CpG_sites_only_LIFR + boxplot_lo_rep_CpG_sites_only_LIF +
  boxplot_lo_rep_CpG_sites_only_MEST + boxplot_lo_rep_CpG_sites_only_NNAT + boxplot_lo_rep_CpG_sites_only_PEG3 +
  boxplot_lo_rep_CpG_sites_only_PEG10 + boxplot_lo_rep_CpG_sites_only_PLAGL1 + boxplot_lo_rep_CpG_sites_only_RTL1 +
  boxplot_lo_rep_CpG_sites_only_SLC2A8 + boxplot_lo_rep_CpG_sites_only_SNRPN + boxplot_lo_rep_CpG_sites_only_SUV39H1 +
  boxplot_lo_rep_CpG_sites_only_TXNIP + boxplot_lo_rep_CpG_sites_only_XIST + guide_area() +
  plot_layout(design = layout_design,
              guides = "collect"
              # axis_titles = "collect_y"
  ) &
  theme(plot.tag = element_text(face = "bold",
                                size = 14,
                                colour = "grey20",
                                hjust = 0.00),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA))
# Save the combined violin plots ----
ggsave("data/analyzed_data/boxplots_20240422_ear_tag_all_samples/20240422_target_combined_boxplot_lo_rep_CpG_sites_only.png",
       plot = combined_boxplots,
       height = 825,
       width = 645,
       dpi = 300,
       unit = "mm")
