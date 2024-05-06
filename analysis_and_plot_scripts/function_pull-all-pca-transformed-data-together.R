# ----
#                    ▓▓░  ▓▓▓  ░▓▓▓
#   -. .-.   .-. .- ▓░ ▓░▓░   ░▓  ░▓ -. .-.   .-. .-
#   ||\|||\ /|||\|| ▓▓▓▓░▓░░▓▓░▓     ||\|||\ /|||\||
#   |/ \|||\|||/ \| ▓░ ▓░▓░ ░▓░▓  ░▓ |/ \|||\|||/ \|
#   ¯   `-¯ `-'   ` ▓░ ▓░ ▓▓▓  ░▓▓▓  ¯   `-¯ `-'   `
#  T H E  A P P L I E D  G E N O M I C S  C E N T R E
#           --- Coding a better harvest ---
# Background ----
# The purpose of this script is to transform the data making:
#   a) the mean % methylation for the whole sample,
#   b) the mean % methylation for each target site of a sample its own column, and
#   c) the % methylation for each CpG site
# This will enable PCAs to be done which will have whole sample methylation, target site metehylation, and each CpG site used as principal components
# This script is **part 2 of 2** for transforming data for PCAs:
#   1) `function_transform-data-for-pca.R`
#   2) `funtion_pull-all-pca-transformed-data-together.R`
# Companion function for a) - Collecting together mean methylation for whole samples ----
pull_all_single_sample_mean_meth_together <- function(input_df) {
  whole_sample_mean_perc_meth <- all_mean_perc_meth_for_pca[0,]
  for (i in unique(input_df$et_index)) {
    mean_meth_single_sample <- get_whole_sample_mean_meth(input_df, i)
    whole_sample_mean_perc_meth <- rbind(whole_sample_mean_perc_meth, mean_meth_single_sample)
  }
  whole_sample_mean_perc_meth <<- whole_sample_mean_perc_meth %>%
    distinct(.keep_all = TRUE)
}
# Companion function for b) - Collecting together mean methylation of all targets for each sample ----
# This is a work in progress as I have 2 major issues that I haven't found a means to take care of in a function:
#   i) How to identify each mean_perc_meth column for which target it is
#   ii) How to preserve all et_index's for every iteration so that I can account for when a sample might be missing data for a target_site