# ----
#                    ▓▓░  ▓▓▓  ░▓▓▓
#   -. .-.   .-. .- ▓░ ▓░▓░   ░▓  ░▓ -. .-.   .-. .-
#   ||\|||\ /|||\|| ▓▓▓▓░▓░░▓▓░▓     ||\|||\ /|||\||
#   |/ \|||\|||/ \| ▓░ ▓░▓░ ░▓░▓  ░▓ |/ \|||\|||/ \|
#   ¯   `-¯ `-'   ` ▓░ ▓░ ▓▓▓  ░▓▓▓  ¯   `-¯ `-'   `
#  T H E  A P P L I E D  G E N O M I C S  C E N T R E
#           --- Coding a better harvest ---
# Background ----
# This short script contains the function for pulling together all the read depth counts data into one data table
# Generate the function ----
pull_all_read_counts_together <- function(input_df) {
  all_read_depth_counts_together <- input_df[0,]  
  for (i in unique(input_df$et_index)) {
    for(j in unique(input_df$target_site)) {
      read_range <- generate_read_depth_counts(input_df, i, j)
      all_read_depth_counts_together <- rbind(all_read_depth_counts_together, read_range)
    }
  }
  # all_read_depth_counts_together_offTarg <<- all_read_depth_counts_together %>%
  all_read_depth_counts_together <<- all_read_depth_counts_together %>%
    select(-c(CpG_site,
              perc_meth,
              meth_C_count,
              unmeth_C_count,
              total_read_counts)) %>%
    distinct(.keep_all = TRUE) #%>%
  # write_csv("data/raw_data/20231130_all_read_depth_counts_together.csv")
}