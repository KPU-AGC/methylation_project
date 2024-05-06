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
generate_read_depth_range <- function(input_df, specific_sample, specific_target) {
  range_counts <<- input_df %>%
    filter(sample == specific_sample & target_site == specific_target) %>%
    select(target_site,
           sample,
           total_count_Cs,
           sample_type) %>%
    mutate(max_count_Cs = max(total_count_Cs)) %>%
    mutate(mean_max_C_count = mean(max_count_Cs)) %>%
    mutate(min_count_Cs = min(total_count_Cs)) %>%
    mutate(mean_min_C_count = mean(min_count_Cs)) %>%
    mutate(range_count_Cs = (max_count_Cs - min_count_Cs)) %>%
    mutate(count_Cs_mean = mean(total_count_Cs)) %>%
    select(target_site,
           sample,
           sample_type,
           mean_max_C_count,
           mean_min_C_count,
           range_count_Cs,
           count_Cs_mean) %>%
    distinct(.keep_all = TRUE)
  # return(range_counts)
}