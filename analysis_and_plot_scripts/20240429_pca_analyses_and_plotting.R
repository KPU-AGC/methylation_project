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
pacman::p_load(factoextra,
               FactoMineR,
               ggcorrplot,
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
# date     2024-04-29                                             
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
source("code/function_transform-data-for-pca.R")
source("code/funtion_pull-all-pca-transformed-data-together.R")
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
# Calculate means of total_counts for each target site ----
pull_all_mean_total_count_together(full_run4_dataset)
# Parse out total_counts that are less than 25% the mean value/representation of CpG sites in the target site ----
pull_all_low_rep_CpG_sites_analyses_together(all_mean_total_count_together)
only_hi_rep_CpG_sites <- all_low_rep_CpG_sites_analyses_together %>%
  filter(low_rep_CpG_site == FALSE) %>%
  filter(total_count > 15)
only_hi_rep_recalc_tot_count_CpG_sites <- pull_all_mean_total_count_together(only_hi_rep_CpG_sites)
# Get the mean total C counts and log transform the data ----
pull_all_mean_perc_meth_and_read_counts_together(only_hi_rep_recalc_tot_count_CpG_sites)
# Re-structure sample list data frame for creation of a covariance matrix ----
sample_list_values_transformed <- sample_list %>%
  mutate(calf_age_values = gsub("newborn",
                                0,
                                age_of_calf_at_collection_rounded_weeks)) %>%
  mutate(calf_age_values = as.numeric(gsub("unknown",
                                           NA,
                                           calf_age_values))) %>%
  mutate(calf_sex_values = gsub("female",
                                1,
                                sex)) %>%
  mutate(calf_sex_values = gsub("male",
                                2,
                                calf_sex_values)) %>%
  mutate(calf_sex_values = as.numeric(gsub("unknown",
                                           NA,
                                           calf_sex_values))) %>%
  mutate(calf_survived_values = gsub("survived_birth",
                                     1,
                              survived_birth)) %>%
  mutate(calf_survived_values = gsub("still_born",
                                     2,
                                     calf_survived_values)) %>%
  mutate(calf_survived_values = gsub("dead_on_arrival",
                                     3,
                                     calf_survived_values)) %>%
  mutate(calf_survived_values = gsub("unknown",
                                     NA,
                                     calf_survived_values)) %>%
  mutate(calf_collection_site_values = gsub("CCF",
                                            1,
                                            collection_site)) %>%
  mutate(calf_collection_site_values = gsub("Bill_Croushore",
                                            2,
                                            calf_collection_site_values)) %>%
  mutate(calf_collection_site_values = gsub("Sicotte",
                                            3,
                                            calf_collection_site_values)) %>%
  mutate(calf_collection_site_values = gsub("Trudeau",
                                            4,
                                            calf_collection_site_values)) %>%
  mutate(calf_collection_site_values = gsub("Westcoast",
                                            5,
                                            calf_collection_site_values)) %>%
  mutate(calf_weight_values = as.numeric(gsub("unknown",
                                              NA,
                                              calf_birth_weight_lbs)))
# Add the metadata to a table for plotting and create a covariance matrix ----
# This quick loop is needed to pull out the mean % methylation for each target while preserving ET index orders for the cases in which there may be missing data for samples
transform_data_for_PCAs <- function(input_df, select_target) {
  selected_column <- input_df %>%
    filter(target_site == select_target) %>%
    select(et_index, mean_perc_meth)
  return(selected_column)
}
transform_data_for_conditions_PCAs <- function(input_df, select_pathology) {
  selected_column <- input_df %>%
    filter(pathology == select_pathology) %>%
    select(et_index, mean_perc_meth)
  return(selected_column)
}
for (i in unique(mean_perc_meth_for_cov_mat$target_site)) {
  assign(paste0("mean_perc_meth_",i),
         transform_data_for_PCAs(mean_perc_meth_for_cov_mat, i))
}
for (i in unique(mean_perc_meth_for_cov_mat$pathology)) {
  assign(paste0("mean_perc_meth_",i),
         transform_data_for_conditions_PCAs(mean_perc_meth_for_cov_mat, i))
}
all_mean_perc_meth_for_pca <- data.frame("et_index" = unique(mean_perc_meth_for_cov_mat$et_index)) %>%
  left_join(mean_perc_meth_DMAP1,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_DMAP1" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_DNMT3A,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_DNMT3A" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_DNMT3B,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_DNMT3B" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_GNAS,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_GNAS" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_H19,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_H19" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_IGF2R,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_IGF2R" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_KCNQ1,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_KCNQ1" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_LIF,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_LIF" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_LIFR,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_LIFR" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_MEST,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_MEST" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_NNAT,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_NNAT" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_PEG3,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_PEG3" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_PEG10,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_PEG10" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_PLAGL1,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_PLAGL1" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_RTL1,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_RTL1" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_SLC2A8,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_SLC2A8" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_SNRPN,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_SNRPN" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_SUV39H1,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_SUV39H1" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_TXNIP,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_TXNIP" = "mean_perc_meth") %>%
  left_join(mean_perc_meth_XIST,
            by = "et_index") %>%
  rename(.,
         "mean_perc_meth_XIST" = "mean_perc_meth") %>%
  left_join(sample_list_values_transformed,
            by = "et_index") %>%
  select(-c(methylation_call_file,
            tube_label,
            calf_number,
            ear_punch_id,
            collection_site,
            sex,
            calf_birth_weight_lbs,
            birth_date,
            date_collection,
            date_received_boviteq,
            age_of_calf_at_collection_rounded_weeks,
            survived_birth
            )) %>%
  write_csv("data/raw_data/20240416_all_mean_perc_meth_for_pca.csv")
colnames(all_mean_perc_meth_for_pca)
column_reorder <- c("et_index",
                    "pathology",
                    "mean_perc_meth_DMAP1",
                    "mean_perc_meth_DNMT3A",
                    "mean_perc_meth_DNMT3B",
                    "mean_perc_meth_GNAS",
                    "mean_perc_meth_H19",
                    "mean_perc_meth_IGF2R",
                    "mean_perc_meth_KCNQ1",
                    "mean_perc_meth_LIF",
                    "mean_perc_meth_LIFR",
                    "mean_perc_meth_MEST",
                    "mean_perc_meth_NNAT",
                    "mean_perc_meth_PEG3",
                    "mean_perc_meth_PEG10",
                    "mean_perc_meth_PLAGL1",
                    "mean_perc_meth_RTL1",
                    "mean_perc_meth_SLC2A8",
                    "mean_perc_meth_SNRPN",
                    "mean_perc_meth_SUV39H1",
                    "mean_perc_meth_TXNIP",
                    "mean_perc_meth_XIST",
                    "calf_age_values",
                    "calf_sex_values",
                    "calf_survived_values",
                    "calf_collection_site_values",
                    "calf_weight_values")
all_mean_perc_meth_for_pca <- all_mean_perc_meth_for_pca[, column_reorder] %>%
  # mutate_if(is.character, as.numeric)
  # replace(.,
  #         is.na(.),
  #         0.1) %>%
  mutate_at(c("mean_perc_meth_DMAP1",
              "mean_perc_meth_DNMT3A",
              "mean_perc_meth_DNMT3B",
              "mean_perc_meth_GNAS",
              "mean_perc_meth_H19",
              "mean_perc_meth_IGF2R",
              "mean_perc_meth_KCNQ1",
              "mean_perc_meth_LIF",
              "mean_perc_meth_LIFR",
              "mean_perc_meth_MEST",
              "mean_perc_meth_NNAT",
              "mean_perc_meth_PEG3",
              "mean_perc_meth_PEG10",
              "mean_perc_meth_PLAGL1",
              "mean_perc_meth_RTL1",
              "mean_perc_meth_SLC2A8",
              "mean_perc_meth_SNRPN",
              "mean_perc_meth_SUV39H1",
              "mean_perc_meth_TXNIP",
              "mean_perc_meth_XIST",
              "calf_age_values",
              "calf_sex_values",
              "calf_survived_values",
              "calf_collection_site_values",
              "calf_weight_values"),
            as.numeric)
# I also want to create a data frame that has just the mean methylation for a whole sample ----
# 2024-04-30 - Okay, so these are something that are really a work in progress
# I have the functions for whole samples down and work well as separate functions
# But each target_site and CpG_site are still something that I need to work out
# As can be seen above, it is still very ungainly to get a data frame for each target site
pull_all_single_sample_mean_meth_together(all_mean_perc_meth_for_pca)
# Create a correlation matrix ----
all_mean_perc_meth_for_pca_cov_mat <- PCA(all_mean_perc_meth_for_pca[,3:27],
                                          scale.unit = TRUE,
                                          graph = FALSE)
whole_sample_mean_perc_meth_cov_mat <- PCA(whole_sample_mean_perc_meth[,23:28],
                                           scale.unit = TRUE,
                                           graph = FALSE)
# Extract the results from the correlation matrix ----
all_mean_perc_meth_for_pca_var <- get_pca_var(all_mean_perc_meth_for_pca_cov_mat)
whole_sample_mean_perc_meth_var <- get_pca_var(whole_sample_mean_perc_meth_cov_mat)
whole_sample_mean_perc_meth_pca_results <- get_pca(whole_sample_mean_perc_meth_cov_mat)
whole_sample_mean_perc_meth_ind <- get_pca_ind(whole_sample_mean_perc_meth_cov_mat)
# Get Eigenvalues and variances from correlation matrix ----
eigenvalues_ts <- get_eigenvalue(all_mean_perc_meth_for_pca_cov_mat)
eigenvalues_ts
eigenvalues_ws <- get_eigenvalue(whole_sample_mean_perc_meth_cov_mat)
eigenvalues_ws
# Generate Scree plot to get an idea of what components to use and how many ----
all_mean_perc_meth_for_pca_scree <- fviz_eig(all_mean_perc_meth_for_pca_cov_mat,
                                             add_labels = TRUE,
                                             ncp = 25,
                                             ylim = c(0,30),
                                             ggtheme = theme_classic())
all_mean_perc_meth_for_pca_scree
whole_sample_mean_perc_meth_scree <- fviz_eig(whole_sample_mean_perc_meth_cov_mat,
                                             add_labels = TRUE,
                                             ncp = 6,
                                             ylim = c(0,30),
                                             ggtheme = theme_classic())
whole_sample_mean_perc_meth_scree
# Generate plots on the contributions of elements ----
whole_sample_mean_perc_meth_contrib <- fviz_contrib(whole_sample_mean_perc_meth_cov_mat,
                                                    choice = "var",
                                                    axes = 2, #contributions of the second dimension
                                                    add_labels = TRUE,
                                                    ylim = c(0,30),
                                                    ggtheme = theme_classic())
whole_sample_mean_perc_meth_contrib
# Generate variable correlation plot to get an idea of what components to use and how many ----
colnames(all_mean_perc_meth_for_pca_var$cos2) <- c("PC1","PC2","PC3","PC4","PC5")
rownames(all_mean_perc_meth_for_pca_var$cos2) <- c("DMAP1","DNMT3A","DNMT3B","GNAS","H19","IGF2R","KCNQ1","LIF","LIFR","MEST","NNAT","PEG3","PEG10","PLAGL1","RTL1","SLC2A8","SNRPN","SUV39H1","TXNIP","XIST","Calf Age","Calf Sex","Calf Survived Birth","Farm Site of Collection","Calf Weight")
list(unique(target_list$target_site))
all_mean_perc_meth_for_pca_corr_plot <- ggcorrplot(all_mean_perc_meth_for_pca_var$cos2,
                                                   method = "circle",
                                                   outline.color = "#ffffff",
                                                   legend.title = "Correlation",
                                                   colors = c("#a8721b",
                                                              "#ffffff",
                                                              "#1b62a8")) +
  theme_classic() +
  xlab("Contributing Variables") +
  ylab("Principal Components") +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(face = "bold",
                                    size = 12,
                                    colour = "grey20"),
        axis.text.x = element_text(face = "bold",
                                   hjust = 0.95,
                                   vjust = 0.95,
                                   size = 10,
                                   angle = 45,
                                   colour = "grey30"),
        axis.title.y = element_text(face = "bold",
                                    size = 12,
                                    colour = "grey20"),
        axis.text.y = element_text(face = "bold",
                                   size = 10,
                                   colour = "grey30"),
        legend.title = element_text(face = "bold",
                                    size = 12,
                                    colour = "grey20"),
        legend.text = element_text(face = "bold",
                                   size = 10,
                                   colour = "grey30"),
        plot.background = element_rect(fill = "transparent",
                                       colour = NA),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA),
        legend.background = element_rect(fill = "transparent",
                                         colour = NA),
        legend.box.background = element_rect(fill = "transparent",
                                             colour = NA))
  
all_mean_perc_meth_for_pca_corr_plot
ggsave("data/analyzed_data/pca_plots_20240430_ear_tag_all_samples/20240430_corr_plot_all_mean_meth_filt_data_only.png",
       width = 10,
       height = 4,
       dpi = 300,
       plot = all_mean_perc_meth_for_pca_corr_plot)
# Okay, so looking at these results of what variables are contributing to principal copmonents it looks like PC1 and PC2 have specific stand outs
# PC1 is majorly made up of:
paste0("DNMT3A (Correlation of ",round(all_mean_perc_meth_for_pca_var$cos2[2,1], digits = 3),")")
paste0("GNAS (Correlation of ",round(all_mean_perc_meth_for_pca_var$cos2[4,1], digits = 3),")")
paste0("H19 (Correlation of ",round(all_mean_perc_meth_for_pca_var$cos2[5,1], digits = 3),")")
paste0("KCNQ1 (Correlation of ",round(all_mean_perc_meth_for_pca_var$cos2[7,1], digits = 3),")")
paste0("MEST (Correlation of ",round(all_mean_perc_meth_for_pca_var$cos2[10,1], digits = 3),")")
paste0("NNAT (Correlation of ",round(all_mean_perc_meth_for_pca_var$cos2[11,1], digits = 3),")")
paste0("PEG3 (Correlation of ",round(all_mean_perc_meth_for_pca_var$cos2[12,1], digits = 3),")")
paste0("PLAGL1 (Correlation of ",round(all_mean_perc_meth_for_pca_var$cos2[14,1], digits = 3),")")
paste0("XIST (Correlation of ",round(all_mean_perc_meth_for_pca_var$cos2[20,1], digits = 3),")")
# DNMT3A (Correlation of 0.472)
# GNAS (Correlation of 0.688)
# H19 (Correlation of 0.411)
# KCNQ1 (Correlation of 0.725)
# MEST (Correlation of 0.377)
# NNAT (Correlation of 0.301)
# PEG3 (Correlation of 0.709)
# PLAGL1 (Correlation of 0.667)
# XIST (Correlation of 0.369)
# PC2 is majorly made up of:
paste0("SUV39H1 (Correlation of ",round(all_mean_perc_meth_for_pca_var$cos2[18,2], digits = 3),")")
paste0("XIST (Correlation of ",round(all_mean_perc_meth_for_pca_var$cos2[20,2], digits = 3),")")
paste0("Calf Sex (Correlation of ",round(all_mean_perc_meth_for_pca_var$cos2[22,2], digits = 3),")")
# SUV39H1 (Correlation of 0.431)
# XIST (Correlation of 0.388)
# Calf Sex (Correlation of 0.702)
# So, I'm going to re-analyze the data with just these variables included and see how that lays out the PCA plots ----
select_variables_for_pca <- all_mean_perc_meth_for_pca %>%
  select(c(et_index,
           pathology,
           mean_perc_meth_DNMT3A,
           mean_perc_meth_GNAS,
           mean_perc_meth_H19,
           mean_perc_meth_KCNQ1,
           mean_perc_meth_MEST,
           mean_perc_meth_NNAT,
           mean_perc_meth_PEG3,
           mean_perc_meth_PLAGL1,
           mean_perc_meth_SUV39H1,
           mean_perc_meth_XIST,
           calf_sex_values))
# Now generate the correlation matrix for the select targets ----
select_variables_cov_mat <- PCA(select_variables_for_pca[,3:13],
                                scale.unit = TRUE,
                                graph = FALSE)
whole_sample_mean_perc_meth_corr_plot <- corrplot(whole_sample_mean_perc_meth_var$cos2,
                                                 is.corr = FALSE)
whole_sample_mean_perc_meth_corr_plot
# I included other variables such as age, sex, etc but I am really interested in % methylation which I see is best represented in dimensions 3 and 4 (Dim.3 and Dim.4, respectively)
all_mean_perc_meth_for_pca_corr_plot_bar <- fviz_contrib(all_mean_perc_meth_for_pca_cov_mat,
                                                         choice = "var",
                                                         axes = 1) # This gives results for the dimensions/principal components
all_mean_perc_meth_for_pca_corr_plot_bar
whole_sample_mean_perc_meth_corr_plot_bar <- fviz_contrib(whole_sample_mean_perc_meth_cov_mat,
                                                         choice = "var",
                                                         axes = 1) # This gives results for the dimensions/principal components
whole_sample_mean_perc_meth_corr_plot_bar
# Another way to represent this and get values for the contributions that variables have:
all_mean_perc_meth_for_pca_dim_desc <- dimdesc(all_mean_perc_meth_for_pca_cov_mat,
                                               axes = c(1:5),
                                               proba = 0.05)
all_mean_perc_meth_for_pca_dim_desc$Dim.5
whole_sample_mean_perc_meth_dim_desc <- dimdesc(whole_sample_mean_perc_meth_cov_mat,
                                               axes = c(1:5),
                                               proba = 0.05)
whole_sample_mean_perc_meth_dim_desc$Dim.4
# Generate variable correlation plot to get an idea of components ----
all_mean_perc_meth_for_pca_var_corr <- fviz_pca_var(all_mean_perc_meth_for_pca_cov_mat,
                                                    col.var = "cos2",
                                                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                    repel = TRUE)
all_mean_perc_meth_for_pca_var_corr
whole_sample_mean_perc_meth_var_corr <- fviz_pca_var(whole_sample_mean_perc_meth_cov_mat,
                                                    col.var = "cos2",
                                                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                    repel = TRUE)
whole_sample_mean_perc_meth_var_corr
# Plot some PCAs ----
all_mean_perc_meth_for_pca_cov_mat$pathology <- as.factor(all_mean_perc_meth_for_pca_cov_mat$pathology,
                                                          levels = c("control",
                                                                     "large_calf",
                                                                     "abnormal",
                                                                     "unknown"))
order_of_pathologies <- c("control" = "Control",
                          "large_calf" = "Large Calf",
                          "abnormal" = "Abnormal",
                          "unknown" = "Not Documented")
pathology_colours <- c("control" = "#006699",
                       "large_calf" = "#f2a40d",
                       "abnormal" = "#840029",
                       "unknown" = "grey40")
select_variables_PCAplot <- fviz_pca_biplot(select_variables_cov_mat,
                                            geom.ind = "point",
                                            pointshape = 21,
                                            fill.ind = select_variables_for_pca$pathology,
                                            addEllipses = TRUE,
                                            label = "none",
                                            col.var = factor(c("mean_perc_meth" = "Methylation",
                                                               "mean_perc_meth" = "Methylation",
                                                               "mean_perc_meth" = "Methylation",
                                                               "mean_perc_meth" = "Methylation",
                                                               "mean_perc_meth" = "Methylation",
                                                               "mean_perc_meth" = "Methylation",
                                                               "mean_perc_meth" = "Methylation",
                                                               "mean_perc_meth" = "Methylation",
                                                               "mean_perc_meth" = "Methylation",
                                                               "mean_perc_meth" = "Methylation",
                                                               "values" = "Metadata")),
                                            repel = TRUE,
                                            legend.title = list(fill = "Calf Condition",
                                                                colour = "Contributions")) +
  theme_classic() +
  xlab(paste0("PC1 (",round(select_variables_cov_mat$eig[1,2], digits = 3),"%)")) +
  ylab(paste0("PC2 (",round(select_variables_cov_mat$eig[2,2], digits = 3),"%)")) +
  geom_text(aes(x = 0,
                y = 3.5,
                label = "Calf Sex"),
            size = 3.5,
            colour = "grey30") +
  geom_text(aes(x = -2.4,
                y = 3,
                label = "SUV39H1"),
            size = 3.5,
            colour = "#ea0049") +
  geom_text(aes(x = 3.3,
                y = 2.5,
                label = "XIST"),
            size = 3.5,
            colour = "#ea0049") +
  geom_text(aes(x = -3.2,
                y = 1.6,
                label = "DNMT3A"),
            size = 3.5,
            colour = "#ea0049") +
  geom_text(aes(x = -3.2,
                y = 0.7,
                label = "H19"),
            size = 3.5,
            colour = "#ea0049") +
  geom_text(aes(x = 3,
                y = 1,
                label = "NNAT"),
            size = 3.5,
            colour = "#ea0049") +
  geom_text(aes(x = 3.3,
                y = 0.7,
                label = "MEST"),
            size = 3.5,
            colour = "#ea0049") +
  geom_text(aes(x = 4,
                y = 0.45,
                label = "GNAS"),
            size = 3.5,
            colour = "#ea0049") +
  geom_text(aes(x = 4.4,
                y = -0.15,
                label = "KCNQ1"),
            size = 3.5,
            colour = "#ea0049") +
  geom_text(aes(x = 2.25,
                y = -0.45,
                label = "PEG3"),
            size = 3.5,
            colour = "#ea0049") +
  scale_fill_manual(breaks = names(order_of_pathologies),
                    values = pathology_colours,
                    name = "Calf Condition",
                    labels = order_of_pathologies) +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(face = "bold",
                                    size = 12,
                                    colour = "grey20"),
        axis.text.x = element_text(face = "bold",
                                   size = 10,
                                   colour = "grey30"),
        axis.title.y = element_text(face = "bold",
                                    size = 12,
                                    colour = "grey20"),
        axis.text.y = element_text(face = "bold",
                                   size = 10,
                                   colour = "grey30"),
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
                                             colour = NA)) +
  # ggpubr::fill_palette(pathology_colours) +
  ggpubr::color_palette(c("#2f3b42",
                          "#ea0049"))
select_variables_PCAplot
ggsave("data/analyzed_data/pca_plots_20240430_ear_tag_all_samples/20240430_pca_plot_select_variables_filt_data_only.png",
        width = 10,
        height = 7,
        dpi = 300,
        plot = select_variables_PCAplot)
whole_sample_mean_meth_PCAplot <- fviz_pca_biplot(whole_sample_mean_perc_meth_cov_mat,
                                                  geom.ind = "point",
                                                  pointshape = 21,
                                                  fill.ind = whole_sample_mean_perc_meth$pathology,
                                                  addEllipses = TRUE,
                                                  label = "none",
                                                  col.var = factor(c("values" = "Metadata",
                                                                     "values" = "Metadata",
                                                                     "values" = "Metadata",
                                                                     "values" = "Metadata",
                                                                     "values" = "Metadata",
                                                                     "mean_perc_meth" = "Methylation")),
                                                  repel = TRUE,
                                                  legend.title = list(fill = "Calf Condition",
                                                                      colour = "Contributions")) +
  theme_classic() +
  xlab(paste0("PC1 (",round(whole_sample_mean_perc_meth_cov_mat$eig[1,2], digits = 3),"%)")) +
  ylab(paste0("PC2 (",round(whole_sample_mean_perc_meth_cov_mat$eig[2,2], digits = 3),"%)")) +
  geom_text(aes(x = 1.5,
                y = 2,
                label = "Calf Survived Birth"),
            size = 3.5,
            colour = "grey30") +
  geom_text(aes(x = 2.5,
                y = 1.5,
                label = "Calf Weight"),
            size = 3.5,
            colour = "grey30") +
  geom_text(aes(x = 0.7,
                y = -1.5,
                label = "Calf Sex"),
            size = 3.5,
            colour = "grey30") +
  geom_text(aes(x = -2.85,
                y = 0.5,
                label = "Calf Age"),
            size = 3.5,
            colour = "grey30") +
  geom_text(aes(x = -3.1,
                y = 0.9,
                label = "Farm Site of Collection"),
            size = 3.5,
            colour = "grey30") +
  geom_text(aes(x = 1.4,
                y = -2.3,
                label = "Differential\nMethylation"),
            size = 3.5,
            colour = "#ea0049") +
  scale_fill_manual(breaks = names(order_of_pathologies),
                    values = pathology_colours,
                    name = "Calf Condition",
                    labels = order_of_pathologies) +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(face = "bold",
                                    size = 12,
                                    colour = "grey20"),
        axis.text.x = element_text(face = "bold",
                                   size = 10,
                                   colour = "grey30"),
        axis.title.y = element_text(face = "bold",
                                    size = 12,
                                    colour = "grey20"),
        axis.text.y = element_text(face = "bold",
                                   size = 10,
                                   colour = "grey30"),
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
                                             colour = NA)) +
  # ggpubr::fill_palette(pathology_colours) +
  ggpubr::color_palette(c("#2f3b42",
                          "#ea0049"))
whole_sample_mean_meth_PCAplot
ggsave("data/analyzed_data/pca_plots_20240430_ear_tag_all_samples/20240430_pca_plot_whole_sample_filt_data_only.png",
       width = 10,
       height = 7,
       dpi = 300,
       plot = whole_sample_mean_meth_PCAplot)
