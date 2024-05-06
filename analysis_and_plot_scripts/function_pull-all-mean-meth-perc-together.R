# ----
#                    ▓▓░  ▓▓▓  ░▓▓▓
#   -. .-.   .-. .- ▓░ ▓░▓░   ░▓  ░▓ -. .-.   .-. .-
#   ||\|||\ /|||\|| ▓▓▓▓░▓░░▓▓░▓     ||\|||\ /|||\||
#   |/ \|||\|||/ \| ▓░ ▓░▓░ ░▓░▓  ░▓ |/ \|||\|||/ \|
#   ¯   `-¯ `-'   ` ▓░ ▓░ ▓▓▓  ░▓▓▓  ¯   `-¯ `-'   `
#  T H E  A P P L I E D  G E N O M I C S  C E N T R E
#           --- Coding a better harvest ---
# Background ----
# This short script contains the function for pulling together all the average percentages of methylation data into one data table
# Generate the function ----
pull_all_mean_perc_means_together <- function(input_df) {
  all_mean_perc_meth_together <- input_df[0,]  
  for (i in unique(input_df$et_index)) {
    for (j in unique(input_df$target_site)) {
      mean_perc_meth <- generate_mean_perc_meth(input_df, i, j)
      all_mean_perc_meth_together <- rbind(all_mean_perc_meth_together, mean_perc_meth)
    }
  }
  all_mean_perc_meth_together <<- all_mean_perc_meth_together %>%
    select(-c(CpG_site,
              perc_meth,
              meth_C_count,
              unmeth_C_count)) %>%
    distinct(.keep_all = TRUE) %>%
    write_csv("data/raw_data/20231202_all_mean_perc_meth_together.csv")
}