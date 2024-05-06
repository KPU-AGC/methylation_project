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
# This script is **part 1 of 2** for transforming data for PCAs:
#   1) `function_transform-data-for-pca.R`
#   2) `funtion_pull-all-pca-transformed-data-together.R`
# Function a) - Whole sample mean methylation ----
get_whole_sample_mean_meth <- function(input_df, select_sample) {
  specific_et_index_whole_sample_mean_perc_meth <- input_df %>%
    filter(et_index == select_sample) %>%
    mutate(whole_sample_mean_perc_meth = rowMeans(.[,3:22], na.rm = TRUE))
  return(specific_et_index_whole_sample_mean_perc_meth)
}
# Function b) - Each target site mean methylation ----
get_ea_target_mean_meth <- function(input_df, select_target) {
  selected_column <- input_df %>%
    filter(target_site == select_target) %>%
    select(et_index, mean_perc_meth)
  return(selected_column)
}
# Function c) - Each CpG site mean methylation ----
# I'm still on the fence about implementing this one. It's going to take some serious figuring out to handle this