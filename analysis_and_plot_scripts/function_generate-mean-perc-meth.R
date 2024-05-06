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
generate_mean_perc_meth <- function(input_df, specific_sample, specific_target) {
  perc_meth_mean_for_target <- input_df %>%
    filter(et_index == specific_sample & target_site == specific_target) %>%
    mutate(mean_perc_meth = mean(perc_meth)) %>%
    mutate(min_perc_meth = min(perc_meth)) %>%
    mutate(max_perc_meth = max(perc_meth)) %>%
    mutate(range_perc_meth = max_perc_meth - min_perc_meth) %>%
    distinct(.keep_all = TRUE)
  return(perc_meth_mean_for_target)
}