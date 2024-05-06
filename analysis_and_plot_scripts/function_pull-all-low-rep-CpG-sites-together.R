# ----
#                    ▓▓░  ▓▓▓  ░▓▓▓
#   -. .-.   .-. .- ▓░ ▓░▓░   ░▓  ░▓ -. .-.   .-. .-
#   ||\|||\ /|||\|| ▓▓▓▓░▓░░▓▓░▓     ||\|||\ /|||\||
#   |/ \|||\|||/ \| ▓░ ▓░▓░ ░▓░▓  ░▓ |/ \|||\|||/ \|
#   ¯   `-¯ `-'   ` ▓░ ▓░ ▓▓▓  ░▓▓▓  ¯   `-¯ `-'   `
#  T H E  A P P L I E D  G E N O M I C S  C E N T R E
#           --- Coding a better harvest ---
# Background ----
# This short script contains the function for pulling together all low_rep_CpG analyses
# It pulls directly upon `function_filter-low-rep-CpG-sites.R`
# *This script is **part 4 of 4** for filtering out false CpG sites from methylation call data using the TABS protocol*
#   1) `function_generate-mean-total-count-CpG-sites.R`
#   2) `function_pull-all-mean-total-count-together.R`
#   3) `function_filter-low-rep-CpG-sites.R`
#   4) `function_pull-all-low-rep-CpG-sites-together.R`
# Generate the function ----
pull_all_low_rep_CpG_sites_analyses_together <- function(input_df) {
  all_low_rep_CpG_sites_analyses_together <- input_df[0,]  
  for (i in unique(input_df$target_site)) {
    for (j in unique(input_df$CpG_site)) {
      specific_low_rep_CpG_site <- filter_low_rep_CpG_sites(input_df, i, j)
      all_low_rep_CpG_sites_analyses_together <- rbind(all_low_rep_CpG_sites_analyses_together, specific_low_rep_CpG_site)
    }
  }
  all_low_rep_CpG_sites_analyses_together <<- all_low_rep_CpG_sites_analyses_together %>%
    distinct(.keep_all = TRUE) # %>%
    # write_csv("data/raw_data/20240422_all_low_rep_CpG_sites_analyses_together.csv")
}