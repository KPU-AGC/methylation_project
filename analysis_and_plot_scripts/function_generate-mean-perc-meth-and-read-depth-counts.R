# ----
#                    ▓▓░  ▓▓▓  ░▓▓▓
#   -. .-.   .-. .- ▓░ ▓░▓░   ░▓  ░▓ -. .-.   .-. .-
#   ||\|||\ /|||\|| ▓▓▓▓░▓░░▓▓░▓     ||\|||\ /|||\||
#   |/ \|||\|||/ \| ▓░ ▓░▓░ ░▓░▓  ░▓ |/ \|||\|||/ \|
#   ¯   `-¯ `-'   ` ▓░ ▓░ ▓▓▓  ░▓▓▓  ¯   `-¯ `-'   `
#  T H E  A P P L I E D  G E N O M I C S  C E N T R E
#           --- Coding a better harvest ---
# Background ----
# This short script contains the function for calculating the average percent methylation of our dataframes
# Generate the function ----
generate_mean_perc_meth_and_read_counts <- function(input_df, specific_sample, specific_target) {
  perc_meth_mean_and_read_counts_for_target <- input_df %>%
    filter(et_index == specific_sample & target_site == specific_target) %>%
    mutate(mean_perc_meth = mean(perc_meth)) %>%
    mutate(min_perc_meth = min(perc_meth)) %>%
    mutate(max_perc_meth = max(perc_meth)) %>%
    mutate(range_perc_meth = max_perc_meth - min_perc_meth) %>%
    mutate(max_read_count = max(total_count)) %>%
    mutate(mean_max_read_count = mean(max_read_count)) %>%
    mutate(min_read_count = min(total_count)) %>%
    mutate(mean_min_read_count = mean(min_read_count)) %>%
    mutate(range_read_count = (max_read_count - min_read_count)) %>%
    mutate(mean_total_count = mean(total_count)) %>%
    mutate(log_of_mean_total_count = log10(mean_total_count)) %>%
    distinct(.keep_all = TRUE)
  return(perc_meth_mean_and_read_counts_for_target)
}