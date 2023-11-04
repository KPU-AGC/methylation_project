library(tidyverse)
library(dplyr)
library(ggplot)

setwd("E:/projects/boviteq/results/20231002_eartag-analysis/")

CpG_df <- read.csv("R-friendly-3.csv")

target_sites <- unique(CpG_df$primer)
samples <- unique(CpG_df$sample)

control_samples <- unique((CpG_df %>% filter(phenotype=="CONTROL"))$sample)
large_offspring_samples <- unique((CpG_df %>% filter(large_calf==1))$sample)

by_position_analysis <- function() {
  
  p_value_per_position <- data.frame(
    sample=character(),
    position=character(),
    CpG_num=integer(),
    primer=character(),
    p_val=double())
  
  df_to_test <- CpG_df %>%
    filter(avg_depth_per_site_per_sample > 10) %>%
    filter(pos_perc_of_avg_depth > 0.5)
  
  for (target_sample in large_offspring_samples) {
    for (target_region in target_sites) {
      
      CpG_sites_in_region <- unique((df_to_test %>% filter(primer==target_region))$chr)
      
      for (CpG_site in CpG_sites_in_region) {
        
        if (nrow(df_to_test %>%
                 filter(primer==target_region) %>%
                 filter(chr==CpG_site) %>%
                 filter(sample==target_sample)) == 0) {
          next
        }
        
        sample_meth_at_site <- (df_to_test %>%
          filter(primer==target_region) %>%
          filter(chr==CpG_site) %>%
          filter(sample==target_sample))$meth_perc
        

        
        control_mean_at_site <- (df_to_test %>%
          filter(primer==target_region) %>%
          filter(phenotype=="CONTROL") %>%
          filter(chr==CpG_site))$meth_perc
        
        mean_data <- mean(control_mean_at_site)
        sd_data <- sd(control_mean_at_site)
        
        z_score <- (sample_meth_at_site - mean_data) / sd_data
        p_value <- 2 * (1 - pnorm(abs(z_score)))
        
        p_value_per_position <- p_value_per_position %>% add_row(
          sample=target_sample,
          position=CpG_site,
          CpG_num=unique((df_to_test %>%
            filter(primer==target_region) %>%
            filter(phenotype=="CONTROL") %>%
            filter(chr==CpG_site))$CpG_num)[1],
          primer=target_region,
          p_val=p_value
        )
      }
    }
  }
  p_value_per_position$adjusted_p_val <- p.adjust(p_value_per_position$p_val, method = "hochberg")
  return(p_value_per_position)
}

by_position_analysis()

write.csv(by_position_analysis(), "testtesttesttest.csv")

  for (i in unique(input_df$chr)) {
    df_to_test <- input_df %>%
      filter(chr==i) %>%
      mutate(meth_perc=meth_perc/100) %>%
      filter(read_count > min_cov)
    x <- filter(df_to_test, phenotype=="CONTROL")$meth_perc
    y <- filter(df_to_test, phenotype=="HAS_ABNORMALITY")$meth_perc
    
    if (length(x) > 15 && length(y) > 15) {
      test_result <- wilcox.test(x, y, exact = FALSE)
      #print(paste(i, unique(df_to_test$CpG_num), unique(df_to_test$primer), test_result$p.value))
      p_value_per_position <- p_value_per_position %>% add_row(
        position=i,
        CpG_num=unique(df_to_test$CpG_num)[1],
        primer=unique(df_to_test$primer)[1],
        p_val=test_result$p.value
      )
    }
  }
  p_value_per_position$p_val <- p.adjust(p_value_per_position$p_val, method = "hochberg")
  print(p_value_per_position[order(p_value_per_position$p_val),])
}
