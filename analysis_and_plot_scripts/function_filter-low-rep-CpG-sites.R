# ----
#                    ▓▓░  ▓▓▓  ░▓▓▓
#   -. .-.   .-. .- ▓░ ▓░▓░   ░▓  ░▓ -. .-.   .-. .-
#   ||\|||\ /|||\|| ▓▓▓▓░▓░░▓▓░▓     ||\|||\ /|||\||
#   |/ \|||\|||/ \| ▓░ ▓░▓░ ░▓░▓  ░▓ |/ \|||\|||/ \|
#   ¯   `-¯ `-'   ` ▓░ ▓░ ▓▓▓  ░▓▓▓  ¯   `-¯ `-'   `
#  T H E  A P P L I E D  G E N O M I C S  C E N T R E
#           --- Coding a better harvest ---
# Background ----
# This script contains the function which filters out CpG sites that have less than 25% the value of the mean_total_counts for its target
# *This script is **part 3 of 4** for filtering out false CpG sites from methylation call data using the TABS protocol*
#   1) `function_generate-mean-total-count-CpG-sites.R`
#   2) `function_pull-all-mean-total-count-together.R`
#   3) `function_filter-low-rep-CpG-sites.R`
#   4) `function_pull-all-low-rep-CpG-sites-together.R`
# This script was originally designed to be used in conjunction with `function_generate-mean-total-count-CpG-sites.R` and `function_pull-all-mean-total-count-together.R`
# The function here relies on the combined output from these two functions above
# The output of this function is a dataframe that has added a column `low_rep_CpG_site` marking whether it is `TRUE` or `FALSE` 
# Generate the function ----
filter_low_rep_CpG_sites <- function(input_df, specific_target, specific_CpG) {
  low_rep_CpG_analysis_for_site_and_target <- input_df %>%
    filter(target_site == specific_target & CpG_site == specific_CpG) %>%
    mutate(low_rep_CpG_site = ifelse(total_count < ((mean_total_count)*0.25), TRUE, FALSE)) %>%
    distinct(.keep_all = TRUE)
  return(low_rep_CpG_analysis_for_site_and_target)
}