# ----
#                    ▓▓░  ▓▓▓  ░▓▓▓
#   -. .-.   .-. .- ▓░ ▓░▓░   ░▓  ░▓ -. .-.   .-. .-
#   ||\|||\ /|||\|| ▓▓▓▓░▓░░▓▓░▓     ||\|||\ /|||\||
#   |/ \|||\|||/ \| ▓░ ▓░▓░ ░▓░▓  ░▓ |/ \|||\|||/ \|
#   ¯   `-¯ `-'   ` ▓░ ▓░ ▓▓▓  ░▓▓▓  ¯   `-¯ `-'   `
#  T H E  A P P L I E D  G E N O M I C S  C E N T R E
#           --- Coding a better harvest ---
# Background ----
# The purpose of this script is to be part of an automated and, relatively, unbiased means of filtering CpG sites from our data
# This is meant to to remove the ones that are not informative or just adding noise or are muting differences that may be in our data
# This script compares a given CpG site across all samples and tests for signficance using analysis of variance tests
# If significance is found for a given CpG site it is kept. The rationale for this is that significance indicates that there is a difference at this site across samples which may be a difference between different conditions (i.e., healthy vs LOS)
# This specific script here is a function that makes up the second portion of the filtering and iteratively pulls the significant results
# Hochberg procedure for multiple comparisons correction was chosen here to:
# 1) avoid Type II errors which overly conservative Bonferroni might give,
# 2) address that our hypotheses are independent, i.e., every iteration of testing a CpG site is its own test, which the Bonferroni-Holm does not make this assumption, and
# 3) we do not need the strictness of a stepwise procedure that a Hommel correction would give
# Generate the function ----
iteratively_find_sites <- function(input_df) {
  # Function to perform iteration
  significance_results <- input_df[0,] %>% 
    mutate(p_value = NA) %>%
    mutate(p_value = NA) %>%
    mutate(p_adjusted = NA) %>%
    mutate(significant = NA) %>%
    mutate(method = NA) %>%
    mutate(normal = NA) %>%
    mutate(equal_variances = NA)
  for (i in unique(input_df$CpG_site)) {
    result <- filtered_CpG_sites_by_ANOVA(input_df, i)
    significance_results <- rbind(significance_results, result)
  }
  sig_w_p_adjusted <<- significance_results %>%
  filter(significant == TRUE) %>%
  adjust_pvalue(.,
                p.col = "p_value",
                output.col = "p_adjusted",
                method = "hochberg") %>%
  mutate(padj_significant = ifelse(p_adjusted >= 0.05, FALSE, TRUE)) %>%
  filter(padj_significant == TRUE)
}