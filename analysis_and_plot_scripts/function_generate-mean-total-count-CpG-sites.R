# ----
#                    ▓▓░  ▓▓▓  ░▓▓▓
#   -. .-.   .-. .- ▓░ ▓░▓░   ░▓  ░▓ -. .-.   .-. .-
#   ||\|||\ /|||\|| ▓▓▓▓░▓░░▓▓░▓     ||\|||\ /|||\||
#   |/ \|||\|||/ \| ▓░ ▓░▓░ ░▓░▓  ░▓ |/ \|||\|||/ \|
#   ¯   `-¯ `-'   ` ▓░ ▓░ ▓▓▓  ░▓▓▓  ¯   `-¯ `-'   `
#  T H E  A P P L I E D  G E N O M I C S  C E N T R E
#           --- Coding a better harvest ---
# Background ----
# This short script contains the function for calculating the average of total counts (methylated or not) for each target
# The purpose of this function is to be part of trimming out CpG sites that are under-represented in the data (i.e., CpG sites that are less than a quarter of the mean number of counts for a site)
# *This script is **part 1 of 4** for filtering out false CpG sites from methylation call data using the TABS protocol*
#   1) `function_generate-mean-total-count-CpG-sites.R`
#   2) `function_pull-all-mean-total-count-together.R`
#   3) `function_filter-low-rep-CpG-sites.R`
#   4) `function_pull-all-low-rep-CpG-sites-together.R`
# **This has a critical distinction from `function_generate-read-depth-counts.R`**
#   **that script is for the purpose of getting inferences on read depth but this script is to aid in parsing out false sites**
# This script was originally designed to be used in conjunction with `function_pull-all-mean-total-count-together.R` which pulls on this function for the complete output
# The output from this script and `function_pull-all-mean-total-count-together.R` is meant to then be used with `function_filter-low-rep-CpG-sites.R`
#   `function_filter-low-rep-CpG-sites.R` performs the actual behaviour of filtering out CpG sites that have less than 25% the value of the mean_total_counts for its target
# Generate the function ----
generate_mean_total_count <- function(input_df, specific_target) {
  total_count_mean_for_target <- input_df %>%
    filter(target_site == specific_target) %>%
    mutate(mean_total_count = mean(total_count)) %>%
    mutate(min_total_count = min(total_count)) %>%
    mutate(max_total_count = max(total_count)) %>%
    mutate(range_total_count = max_total_count - min_total_count) %>%
    distinct(.keep_all = TRUE)
  return(total_count_mean_for_target)
}