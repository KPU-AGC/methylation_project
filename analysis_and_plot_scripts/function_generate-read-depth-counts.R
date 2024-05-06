# ----
#                    ▓▓░  ▓▓▓  ░▓▓▓
#   -. .-.   .-. .- ▓░ ▓░▓░   ░▓  ░▓ -. .-.   .-. .-
#   ||\|||\ /|||\|| ▓▓▓▓░▓░░▓▓░▓     ||\|||\ /|||\||
#   |/ \|||\|||/ \| ▓░ ▓░▓░ ░▓░▓  ░▓ |/ \|||\|||/ \|
#   ¯   `-¯ `-'   ` ▓░ ▓░ ▓▓▓  ░▓▓▓  ¯   `-¯ `-'   `
#  T H E  A P P L I E D  G E N O M I C S  C E N T R E
#           --- Coding a better harvest ---
# Background ----
# This short script contains the function for calculating the range of read depth using datasets from all NGS runs
# The minimum is a average of the minimum read depth across CpG sites for a target gene
# The maximum is a average of the maximum read depth across CpG sites for a target gene
# Generate the function ----
generate_read_depth_counts <- function(input_df, specific_sample, specific_target) {
  range_counts <- input_df %>%
    filter(et_index == specific_sample & target_site == specific_target) %>%
    mutate(total_read_counts = meth_C_count + unmeth_C_count) %>%
    mutate(max_read_count = max(total_read_counts)) %>%
    mutate(mean_max_read_count = mean(max_read_count)) %>%
    mutate(min_read_count = min(total_read_counts)) %>%
    mutate(mean_min_read_count = mean(min_read_count)) %>%
    mutate(range_read_count = (max_read_count - min_read_count)) %>%
    mutate(mean_total_count = mean(total_read_counts)) %>%
    mutate(log_of_mean_total_count = log10(mean_total_count)) %>%
    distinct(.keep_all = TRUE)
  return(range_counts)
}