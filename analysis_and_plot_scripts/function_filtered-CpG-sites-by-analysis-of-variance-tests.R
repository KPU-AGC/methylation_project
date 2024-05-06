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
# This specific script here is a function that makes up the component of the script that filters the CpG sites
# Generate the function ----
filtered_CpG_sites_by_ANOVA <- function(input_df, CpG_position) {
  # Function filters CpG site to a given position and performs ANOVA depending
  # on tests of normality and equal variances. Returns a dataframe filtered to
  # the given position with statistical significance and results of tests of
  # normality and equal variances
  # 
  # Parameters:
  #   input_df (dataframe): input dataframe
  #   CpG_position (chr): CpG position for filtering
  #
  # Returns:
  #   (dataframe): filtered dataframe with test results
  test_data <- input_df %>%
    filter(CpG_site == CpG_position) %>%
    group_by(pathology) %>%
    filter(n() >= 2) %>%
    ungroup()
  if (length(unique(test_data$pathology)) <= 2) {
    test_data <- test_data %>% 
      mutate(p_value = NA) %>%
      mutate(p_adjusted = NA) %>%
      mutate(significant = NA) %>%
      mutate(method = NA) %>%
      mutate(normal = NA) %>%
      mutate(equal_variances = NA)
    return(test_data)
  }
  linear_model_CpGs <- lm(perc_meth ~ pathology,
                          data = test_data)
  shapiroTest <- shapiro_test(linear_model_CpGs$residuals)
  levenesTest <- test_data %>%
    mutate(pathology = factor(pathology)) %>%
    levene_test(perc_meth ~ pathology)
  if (shapiroTest$p.value >= 0.05 & levenesTest$p >= 0.05) {
    normality <- TRUE
    equal_variances <- TRUE
    test_results <- test_data %>%
      anova_test(perc_meth ~ pathology)
    method = "ANOVA"
  } else if (shapiroTest$p.value >= 0.05 & levenesTest$p <= 0.05) {
    normality <- TRUE
    equal_variances <- FALSE
    test_results <- test_data %>%
      welch_anova_test(perc_meth ~ pathology)
    method = "Welchs_ANOVA"
  } else if (shapiroTest$p.value <= 0.05 & levenesTest$p >= 0.05) {
    normality <- FALSE
    equal_variances <- TRUE
    test_results <- test_data %>%
      kruskal_test(perc_meth ~ pathology)
    method = "Kruskal-Wallis"
  } else {
    normality <- FALSE
    equal_variances <- FALSE
    test_results <- test_data %>%
      kruskal_test(perc_meth ~ pathology)
    method = "Kruskal-Wallis"
  }
  test_data <- test_data %>% 
    mutate(p_value = test_results$p) %>%
    mutate(p_adjusted = NA) %>%
    mutate(significant = ifelse(test_results$p >= 0.05, FALSE, TRUE)) %>%
    mutate(method = method) %>%
    mutate(normal = normality) %>%
    mutate(equal_variances = equal_variances)
  return(test_data)
}