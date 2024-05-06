# ----
#                    ▓▓░  ▓▓▓  ░▓▓▓
#   -. .-.   .-. .- ▓░ ▓░▓░   ░▓  ░▓ -. .-.   .-. .-
#   ||\|||\ /|||\|| ▓▓▓▓░▓░░▓▓░▓     ||\|||\ /|||\||
#   |/ \|||\|||/ \| ▓░ ▓░▓░ ░▓░▓  ░▓ |/ \|||\|||/ \|
#   ¯   `-¯ `-'   ` ▓░ ▓░ ▓▓▓  ░▓▓▓  ¯   `-¯ `-'   `
#  T H E  A P P L I E D  G E N O M I C S  C E N T R E
#           --- Coding a better harvest ---
# Background ----
# The purpose of this script is to test for significant differences in read depth between our different pathologies in the ear tag data
# It's intended for use with data that's been transformed to have average percent methylation calculated for a given target in each sample
# This is also intended for use in making target-to-target comparisons for all ear tag samples grouped into their pathologies
# Generate the function ----
test_mean_total_read_count_by_ANOVA <- function(input_df) {
  test_data <- input_df %>%
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
  linear_model_read_count <- lm(mean_total_count ~ pathology,
                               data = select_target_to_plot)
  shapiroTest <- shapiro_test(linear_model_read_count$residuals)
  LevenesTest <- select_target_to_plot %>%
    mutate(pathology = factor(pathology)) %>%
    levene_test(mean_total_count ~ pathology)
  if (shapiroTest$p.value >= 0.05 & LevenesTest$p >= 0.05) {
    normality <- TRUE
    equal_variances <- TRUE
    test_results <- test_data %>%
      anova_test(mean_total_count ~ pathology)
    method = "ANOVA"
  } else if (shapiroTest$p.value >= 0.05 & LevenesTest$p <= 0.05) {
    normality <- TRUE
    equal_variances <- FALSE
    test_results <- test_data %>%
      welch_anova_test(mean_total_count ~ pathology)
    method = "Welchs_ANOVA"
  } else if (shapiroTest$p.value <= 0.05 & LevenesTest$p >= 0.05) {
    normality <- FALSE
    equal_variances <- TRUE
    test_results <- test_data %>%
      kruskal_test(mean_total_count ~ pathology)
    method = "Kruskal-Wallis"
  } else {
    normality <- FALSE
    equal_variances <- F
    test_results <- test_data %>%
      kruskal_test(mean_total_count ~ pathology)
    method = "Kruskal-Wallis"
  }
  test_data <<- test_data %>% 
    mutate(p_value = test_results$p) %>%
    mutate(p_adjusted = NA) %>%
    mutate(significant = ifelse(test_results$p >= 0.05, FALSE, TRUE)) %>%
    mutate(method = method) %>%
    mutate(normal = normality) %>%
    mutate(equal_variances = equal_variances)
}