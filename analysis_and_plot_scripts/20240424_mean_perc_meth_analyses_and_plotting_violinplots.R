# ----
#                    ▓▓░  ▓▓▓  ░▓▓▓
#   -. .-.   .-. .- ▓░ ▓░▓░   ░▓  ░▓ -. .-.   .-. .-
#   ||\|||\ /|||\|| ▓▓▓▓░▓░░▓▓░▓     ||\|||\ /|||\||
#   |/ \|||\|||/ \| ▓░ ▓░▓░ ░▓░▓  ░▓ |/ \|||\|||/ \|
#   ¯   `-¯ `-'   ` ▓░ ▓░ ▓▓▓  ░▓▓▓  ¯   `-¯ `-'   `
#  T H E  A P P L I E D  G E N O M I C S  C E N T R E
#           --- Coding a better harvest ---
# Background ----
# This script is for plotting out Boviteq NGS run 4 assessing ear tag samples
# This is yet again another updated approach with polished code and plots
# I know have a script that filters  
# The purpose of this script is to plot the mean percent methylation of individual samples for a given target where samples are grouped by their pathology/condition
# Loading packages ----
library(pacman)
pacman::p_load(patchwork,
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
# date     2024-04-27                                             
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
# Source functions ----
source("~/Documents/general_code/getOutliers_function.R")
source("code/function_generate-mean-total-count-CpG-sites.R")
source("code/function_pull-all-mean-total-count-together.R")
source("code/function_filter-low-rep-CpG-sites.R")
source("code/function_pull-all-low-rep-CpG-sites-together.R")
source("code/function_generate-mean-perc-meth-and-read-depth-counts.R")
source("code/function_pull-all-mean-meth-perc-and-read-counts-together.R")
source("code/function_test-pathologies-by-analysis-of-variance.R")
# Loading Metadata ----
sample_list <- read_csv("data/analyzed_data/metadata_formatted_for_parsing_single_pathologies.csv")
target_list <- read_csv("data/analyzed_data/bsx_target_sites_trimmed_target_name.csv") %>%
  rename("target_site" = "primer")
# Loading Data ----
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
# Calculate means of total_counts for each target site ----
pull_all_mean_total_count_together(full_run4_dataset)
# Parse out total_counts that are less than 25% the mean value/representation of CpG sites in the target site ----
pull_all_low_rep_CpG_sites_analyses_together(all_mean_total_count_together)
only_hi_rep_CpG_sites <- all_low_rep_CpG_sites_analyses_together %>%
  filter(low_rep_CpG_site == FALSE) %>%
  filter(total_count > 15) 
# Get the mean total C counts and log transform the data ----
pull_all_mean_perc_meth_and_read_counts_together(only_hi_rep_CpG_sites)
# Add the metadata to a table for plotting ----
only_hi_rep_recalc_tot_count_CpG_sites_for_plots <- left_join(all_mean_perc_meth_and_read_counts_together,
                                                              sample_list,
                                                              by = "et_index") #%>%
# filter(!pathology %in% c("abnormal","unknown"))
length(unique(all_mean_perc_meth_and_read_counts_together$et_index))
# [1] 64
length(unique(only_hi_rep_recalc_tot_count_CpG_sites_for_plots$et_index))
# [1] 64
# Universal odds and ends for the plots ----
order_of_sexes <- c("female" = "Female",
                    "male" = "Male",
                    "unknown" = "Not Recorded")
sex_colours <- c("female" = "#209e2c",
                 "male" = "#841ab5",
                 "unknown" = "grey20")
order_of_pathologies <- c("control" = "Control",
                          "large_calf" = "Large Calf",
                          "abnormal" = "Abnormal",
                          "unknown" = "Not Documented")
pathology_colours <- c("control" = "#006699",
                       "large_calf" = "#f2a40d",
                       "abnormal" = "#840029",
                       "unknown" = "grey40")
# Plotting the percent methylation and read depth range ----
for (i in unique(only_hi_rep_recalc_tot_count_CpG_sites_for_plots$target_site)) {
  # Select the specific target and determine outliers ----
  select_target_to_plot <- only_hi_rep_recalc_tot_count_CpG_sites_for_plots %>%
    filter(target_site == i) %>%
    # mutate(outlier = et_index) %>%
    mutate(outlier = ifelse(is_outlier(mean_perc_meth),
                            et_index,
                            as.numeric(NA))) %>%
    mutate(is_Kune = ifelse(et_index == "LLN-Kune",
                            et_index,
                            as.numeric(NA))) %>%
    mutate(mean_perc_meth_wo_KUNE = ifelse(et_index == "LLN-Kune",
                                           as.numeric(NA),
                                           mean_perc_meth))
  # select_target_to_plot$outlier[which(is.na(select_target_to_plot$perc_meth))] <- as.numeric(NA)
  # Perform the stats ----
  test_mean_perc_meth_by_ANOVA(select_target_to_plot)
  # Label the imprinting statuses ----
  imprinted_status <- target_list %>%
    filter(target_site == i) %>%
    mutate_if(is.character, ~replace_na(.,""))
  target_chromosome_number <- target_list %>%
    filter(target_site == i) %>%
    mutate_if(is.numeric, ~replace_na(.,""))
  # Plot the percent methylation ----
  violinplot <- ggplot(data = select_target_to_plot,
                       aes(x = factor(pathology,
                                      levels = c("control",
                                                 "large_calf",
                                                 "abnormal",
                                                 "unknown")),
                           y = mean_perc_meth_wo_KUNE,
                           fill = pathology)) +
    scale_y_continuous(name = "Percent Methylation (%)",
                       limits = c(0,110),
                       breaks = c(0,25, 50, 75, 100)) +
    # sec.axis = sec_axis(~ (./yaxis_range_slope)-yaxis_range_intercept,
    #                     name = "Average Range of Read Depth",
    #                     breaks = seq(0,6,by=1),
    #                     labels = c("1","10","100","1,000","10,000","100,000","1,000,000"))) +
    scale_x_discrete(limits = c("control",
                                "large_calf",
                                "abnormal",
                                "unknown"),
                     breaks = c("control",
                                "large_calf",
                                "abnormal",
                                "unknown"),
                     labels = order_of_pathologies) +
  # geom_linerange(aes(x = pathology,
  #                    ymin = yaxis_range_slope*(log10(mean_min_read_count))+yaxis_range_intercept,
  #                    ymax = yaxis_range_slope*(log10(mean_max_read_count))+yaxis_range_intercept),
  #                linewidth = 12,
  #                alpha = 0.6,
  #                colour = "#16a864") +
    geom_violin(colour = "#2f3b42",
                width = 0.8,
                trim = TRUE,
                scale = "width",
                lwd = 0.2,
                show.legend = TRUE,
                position = position_dodge(1)) +
    stat_summary(fun.data = "mean_cl_boot",
                 geom = "pointrange",
                 colour = "#2f3b42",
                 size = 0.5,
                 show.legend = FALSE) +
    geom_jitter(aes(x = factor(pathology,
                               levels = c("control",
                                          "large_calf",
                                          "abnormal",
                                          "unknown")),
                    y = mean_perc_meth,
                    colour = sex),
                position = position_jitter(0.1),
                size = 0.3,
                show.legend = TRUE) +
    scale_colour_manual(breaks = names(order_of_sexes),
                        values = sex_colours,
                        name = "Calf Sex",
                        labels = order_of_sexes) +
    scale_fill_manual(breaks = names(order_of_pathologies),
                      values = pathology_colours,
                      name = "Calf Condition",
                      labels = order_of_pathologies) +
    geom_segment(aes(x = 1,
                     y = max(mean_perc_meth) + 6,
                     xend = 4,
                     yend = max(mean_perc_meth) + 6),
                 colour = "grey20") +
    geom_segment(aes(x = 1,
                     y = max(mean_perc_meth) + 2,
                     xend = 1,
                     yend = max(mean_perc_meth) + 6),
                 colour = "grey20") +
    geom_segment(aes(x = 4,
                   y = max(mean_perc_meth) + 2,
                   xend = 4,
                   yend = max(mean_perc_meth) + 6),
               colour = "grey20") +
    geom_text(aes(x = 2.5,
                  y = max(mean_perc_meth) + 9,
                  label = paste0("p = ",test_data$p_value)),
              size = 2.5,
              colour = "grey30") +
    theme_classic() +
    xlab("Calf Condition") +
    ggtitle(paste0(i,imprinted_status$imprinted_asterisk)) +
    # geom_text(aes(x = factor(pathology,
    #                          levels = c("control",
    #                                     "large_calf",
    #                                     "abnormal",
    #                                     "unknown")),
    #               y = mean_perc_meth_wo_KUNE,
    #               label = outlier),
    #           na.rm = TRUE,
    #           hjust = -0.1,
    #           size = 3,
    #           position = position_nudge(x = 0.1),
    #           colour = "#ea0049") +
    geom_text(aes(x = factor(pathology,
                             levels = c("control",
                                        "large_calf",
                                        "abnormal",
                                        "unknown")),
                  y = mean_perc_meth,
                  label = is_Kune),
              na.rm = TRUE,
              hjust = -0.1,
              size = 3,
              position = position_nudge(x = 0.1),
              colour = "#ea0049") +
    theme(plot.title = element_text(hjust = 0.5,
                                    face = "bold",
                                    size = 18,
                                    colour = "grey20"),
          axis.title.x = element_text(face = "bold",
                                      size = 12,
                                      colour = "grey20"),
          axis.text.x = element_text(face = "bold",
                                     size = 10,
                                     angle = 45,
                                     hjust = 0.95,
                                     vjust = 0.95,
                                     colour = "grey30"),
          axis.title.y = element_text(face = "bold",
                                      size = 12,
                                      colour = "grey20"),
          axis.text.y = element_text(face = "bold",
                                     size = 10,
                                     colour = "grey30"),
          # axis.title.y.right = element_text(face = "bold",
          #                                   size = 12,
          #                                   vjust = 1,
          #                                   colour = "grey20"),
          # axis.text.y.right = element_text(face = "bold",
          #                                  size = 8,
          #                                  colour = "grey30"),
          legend.title = element_text(face = "bold",
                                      size = 12,
                                      colour = "grey20"),
          legend.text = element_text(face = "bold",
                                     size = 12,
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
  export_plot <- violinplot +
    plot_annotation(caption = paste0("Notes:\n• ",imprinted_status$imprinted_asterisk,"Target is ",imprinted_status$which_parent_imprinted," imprinted\n• Target is on Chromosome ",target_chromosome_number$chromosome_no,"\n• LLN-Kune is labelled if present in Control samples\n• Outliers from the 0.25 to 0.75 quantiles are labelled\n• By ",test_data$method,": p-value = ",test_data$p_value)) &
    theme(plot.caption = element_text(hjust = 0.00))
  # Save the plot ----
  ggsave(paste0("data/analyzed_data/violinplots_20240422_ear_tag_all_samples/20240422_target_",i,"_violinplot_filt_mean_meth_perc_only.png"),
         width = 8,
         height = 7,
         dpi = 300,
         plot = export_plot)
  # assign(paste0("violinplot_filt_mean_meth_perc_only_",i),
  #        violinplot)
}
# Violinplots for poster figures ----
violinplot_filt_mean_meth_perc_only_NNAT
# layout_design <- "
# AABBCCDDE
# FFGGHHII#
# JJKKLLMM#
# NNOOPPQQ#
# RRSSTTUU#
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
combined_violinplots <- violinplot_filt_mean_meth_perc_only_DMAP1 + violinplot_filt_mean_meth_perc_only_DNMT3A + violinplot_filt_mean_meth_perc_only_DNMT3B +
  violinplot_filt_mean_meth_perc_only_GNAS + violinplot_filt_mean_meth_perc_only_H19 + violinplot_filt_mean_meth_perc_only_IGF2R +
  violinplot_filt_mean_meth_perc_only_KCNQ1 + violinplot_filt_mean_meth_perc_only_LIFR + violinplot_filt_mean_meth_perc_only_LIF +
  violinplot_filt_mean_meth_perc_only_MEST + violinplot_filt_mean_meth_perc_only_NNAT + violinplot_filt_mean_meth_perc_only_PEG3 +
  violinplot_filt_mean_meth_perc_only_PEG10 + violinplot_filt_mean_meth_perc_only_PLAGL1 + violinplot_filt_mean_meth_perc_only_RTL1 +
  violinplot_filt_mean_meth_perc_only_SLC2A8 + violinplot_filt_mean_meth_perc_only_SNRPN + violinplot_filt_mean_meth_perc_only_SUV39H1 +
  violinplot_filt_mean_meth_perc_only_TXNIP + violinplot_filt_mean_meth_perc_only_XIST +  guide_area() +
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
# violinplots_for_poster_fig <- violinplot_filt_data_DMAP1 + violinplot_filt_data_DNMT3A + violinplot_filt_data_DNMT3B + violinplot_filt_data_KCNQ1 + guide_area() +
#   violinplot_filt_data_LIF + violinplot_filt_data_NNAT + violinplot_filt_data_RTL1 + violinplot_filt_data_TXNIP +
#   plot_layout(design = layout_design,
#               guides = "collect") &
#   theme(plot.tag = element_text(face = "bold",
#                                 size = 14,
#                                 colour = "grey20",
#                                 hjust = 0.00),
#         plot.background = element_rect(fill = "transparent",
#                                        colour = NA))
# violinplots_for_poster_fig
# ggsave("~/Documents/conferences_events/PAG_2024/posters/assets_resources_figures/20240108_ear_tag_samples_all_pathologies_for_poster.png",
#        plot = violinplots_for_poster_fig,
#        height = 200,
#        width = 400.603,
#        unit = "mm",
#        dpi = 300)
# ggsave("~/Documents/conferences_events/PAG_2024/posters/assets_resources_figures/20240108_ear_tag_samples_all_pathologies_for_poster.pdf",
#        plot = violinplots_for_poster_fig,
#        height = 200,
#        width = 400.603,
#        unit = "mm")
ggsave("data/analyzed_data/violinplots_20240422_ear_tag_all_samples/20240422_combined_violinplots_filt_mean_meth_perc_only.png",
       plot = combined_violinplots,
       height = 850,
       width = 400,
       dpi = 300,
       unit = "mm")
# # Violinplots for direct comparisons ----
# for (i in unique(list_of_targets_w_sig_p$target_site)) {
#   violinplots_for_comparisons <- (paste0("violinplot_",i) | paste0("violinplot_filt_data",i))
#   ggsave(paste0("data/analyzed_data/violinplots_20231211_ear_tag_all_samples_side_by_side_filtered_data_plots/20231211_target_",i,"_side_by_side_comparison.png"),
#          width = 8,
#          height = 7,
#          plot = violinplots_for_comparisons)
# }
