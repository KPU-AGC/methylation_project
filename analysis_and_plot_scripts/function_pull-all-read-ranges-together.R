# ----
#                    ▓▓░  ▓▓▓  ░▓▓▓
#   -. .-.   .-. .- ▓░ ▓░▓░   ░▓  ░▓ -. .-.   .-. .-
#   ||\|||\ /|||\|| ▓▓▓▓░▓░░▓▓░▓     ||\|||\ /|||\||
#   |/ \|||\|||/ \| ▓░ ▓░▓░ ░▓░▓  ░▓ |/ \|||\|||/ \|
#   ¯   `-¯ `-'   ` ▓░ ▓░ ▓▓▓  ░▓▓▓  ¯   `-¯ `-'   `
#  T H E  A P P L I E D  G E N O M I C S  C E N T R E
#           --- Coding a better harvest ---
# Background ----
# This short script contains the function for pulling together all the read depth ranges data into one data table
# Generate the function ----
pull_sites_together <- function(input_df) {
  all_read_depth_ranges_together <- input_df[0,]  
  for (i in unique(input_df$sample)) {
    for(j in unique(input_df$target_site)) {
      read_range <- generate_read_depth_range(input_df, i, j)
      all_read_depth_ranges_together <- rbind(all_read_depth_ranges_together, read_range)
    }
  }
  return(all_read_depth_ranges_together)
}