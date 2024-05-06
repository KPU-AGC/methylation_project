# ----
#                    ▓▓░  ▓▓▓  ░▓▓▓
#   -. .-.   .-. .- ▓░ ▓░▓░   ░▓  ░▓ -. .-.   .-. .-
#   ||\|||\ /|||\|| ▓▓▓▓░▓░░▓▓░▓     ||\|||\ /|||\||
#   |/ \|||\|||/ \| ▓░ ▓░▓░ ░▓░▓  ░▓ |/ \|||\|||/ \|
#   ¯   `-¯ `-'   ` ▓░ ▓░ ▓▓▓  ░▓▓▓  ¯   `-¯ `-'   `
#  T H E  A P P L I E D  G E N O M I C S  C E N T R E
#           --- Coding a better harvest ---
# Background ----
# This script is part of calculating average percent methylation for targets across all samples from all 4 NGS runs
# It works in conjunction with the function pull_sites_together(), which that in turn pulls on the function generate_read_depth_range()
# The purpose of pulling this script below into being its own function is due to the fact that it's quite cumbersome to look at on its own due to all the regex replacements
# Generate the function ----
format_and_generage_read_depth_ranges <- function(list_of_datasets) {
  for (i in 1:length(list_of_datasets)) {
    df_target_site_transformed <- list_of_datasets[[i]] %>%
      select(arb_count,
             target_site,
             sample,
             meth_C_count,
             unmeth_C_count)
    target_site_list_for_replacement <- stri_replace_all_regex(as.list(df_target_site_transformed$target_site),
                                                               pattern = c("NC_037330.1:101[0-9]{6}",
                                                                           "NC_037338.1:74[0-9]{6}",
                                                                           "NC_037340.1:62[0-9]{6}",
                                                                           "NC_037340.1:57[0-9]{6}",
                                                                           "NC_037356.1:49[0-9]{6}",
                                                                           "NC_037336.1:96[0-9]{6}",
                                                                           "NC_037356.1:48[0-9]{6}",
                                                                           "NC_037344.1:69[0-9]{6}",
                                                                           "NC_037347.1:35[0-9]{6}",
                                                                           "NC_037331.1:94[0-9]{6}",
                                                                           "NC_037340.1:66[0-9]{6}",
                                                                           "NC_037331.1:12[0-9]{6}",
                                                                           "NC_037345.1:64[0-9]{6}",
                                                                           "NC_037336.1:81[0-9]{6}",
                                                                           "NC_037348.1:65[0-9]{6}",
                                                                           "NC_037338.1:98[0-9]{6}",
                                                                           "NC_037348.1:19[0-9]{5}",
                                                                           "NC_037357.1:86[0-9]{6}",
                                                                           "NC_037330.1:21[0-9]{6}",
                                                                           "NC_037357.1:77[0-9]{6}"),
                                                               replacement = c("DMAP1",
                                                                               "DNMT3A",
                                                                               "DNMT3B",
                                                                               "GNAS",
                                                                               "H19",
                                                                               "IGF2R",
                                                                               "KCNQ1",
                                                                               "LIF",
                                                                               "LIFR",
                                                                               "MEST",
                                                                               "NNAT",
                                                                               "PEG10",
                                                                               "PEG3",
                                                                               "PLAGL1",
                                                                               "RTL1",
                                                                               "SLC2A8",
                                                                               "SNRPN",
                                                                               "SUV39H1",
                                                                               "TXNIP",
                                                                               "XIST"),
                                                               vectorize_all = FALSE) %>%
      as.data.frame() %>%
      rename("target_site" = ".") %>%
      mutate(arb_count = seq(1,length(target_site)))
    df_target_site_transformed <- df_target_site_transformed %>%
      select(arb_count) %>%
      left_join(.,
                target_site_list_for_replacement,
                by = "arb_count")
    df_minus_target_site <- list_of_datasets[[i]] %>%
      select(!target_site) %>%
      unite("CpG_loci",
            CpG_loci,
            end,
            sep = "-",
            remove = TRUE)
    df_together_w_transformed_target_site <- full_join(df_target_site_transformed,
                                                       df_minus_target_site,
                                                       by = "arb_count",
                                                       relationship = "many-to-many") %>%
      select(!arb_count) %>%
      mutate(total_count_Cs = (meth_C_count + unmeth_C_count)) %>%
      filter(total_count_Cs >= 15) %>%
      left_join(.,
                sample_list,
                by = "sample",
                relationship = "many-to-many") %>%
      select(!Index) %>%
      filter(!grepl("NC_",
                    target_site)) %>%
      filter(!grepl("NW_",
                    target_site)) %>%
      write_csv(paste0("data/raw_data/20231122_run",i,"_together_for_read_range_counts.csv"))
    assign(paste0("run",i,"_together_for_read_range_counts"),
           df_together_w_transformed_target_site)
  }
  # Calculate the read ranges for target sites
  list_of_runs_all_prepped_for_range_calc <- list(run1_together_for_read_range_counts,
                                                  run2_together_for_read_range_counts,
                                                  run3_together_for_read_range_counts,
                                                  run4_together_for_read_range_counts)
  for (i in 1:length(list_of_runs_all_prepped_for_range_calc)) {
    df_read_depth_range <- pull_sites_together(list_of_runs_all_prepped_for_range_calc[[i]])
    assign(paste0("run",i,"_read_depth_range"),
           df_read_depth_range)
  }
  all_runs_together <<- bind_rows(run1_read_depth_range,
                                 run2_read_depth_range,
                                 run3_read_depth_range,
                                 run4_read_depth_range) %>%
    drop_na(sample_type)
}