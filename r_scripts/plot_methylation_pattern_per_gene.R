library(tidyverse)
library(ggplot2)

setwd("E:/projects/boviteq/results/20231002_eartag-analysis/")
CpG_df <- read.csv("methylation_per_sample.csv")

target_sites <- unique(CpG_df$primer)
samples <- unique(CpG_df$sample)

plot_by_gene <- function(
    input_df, 
    input_gene,
    specific_CpG = NULL,
    specific_pos = NULL,
    min_cov = 100,
    show_boxplot = TRUE,
    show_points = TRUE,
    target_sample = NULL,
    plot_chr = FALSE) {
  # Function generates a boxplot of CpG methylation across a given gene. 
  # 
  # Parameters:
  #   input_df (dataframe): input long-format dataframe of methylation status
  #   input_gene (string): input gene for filtering, should be the primer name -- e.g., DNMT3A-BP1X
  #   specific_CpG (int): filter to specific numbered CpG site
  #   specific_pos (float): filter to specific numbered CpG in genomic coordinate form
  #   min_cov (int): filter using minimum coverage
  #   show_boxplot (bool): show boxplots
  #   show_points (bool): show points
  #   target_sample (string): in addition to the control and abnormal boxplots/points, also highlight a given sample
  #   plot_chr (bool): instead of CpG numbering, use genomic coordinates

  jitter_width = 0.2
  jitter_seed = 1234
  
  plotting_column <- ifelse(plot_chr, "chr", "CpG_num")
  
  # filter dataframe according to target gene, minimum coverage, abnormality
  plot_df <- input_df %>%
    filter(primer==input_gene) %>%
    filter(avg_depth_per_site_per_sample > min_cov) %>%
    filter(pos_perc_of_avg_depth > 0.5) %>%
    filter( ((phenotype == "HAS_ABNORMALITY") & (large_calf == 1)) | (phenotype == "CONTROL")) %>%
    mutate(plotting_x = as.factor(get(plotting_column)))

  # deal with plotting target_sample separately
  if ( is.null(target_sample) ) {
    plot_df <- plot_df %>% mutate(visualization_factor=as.factor(ifelse(phenotype=="CONTROL", "CONTROL", "HAS_ABNORMALITY")))
  } else {
    plot_df <- plot_df %>% mutate(visualization_factor=as.factor(if_else(sample==target_sample, target_sample, ifelse(phenotype=="CONTROL", "CONTROL", "HAS_ABNORMALITY"))))
  }
  
  # deal with plotting specific CpG number or genomic coordinate
  if (!is.null(specific_CpG)) {
    plot_df <- plot_df %>% filter(CpG_num==specific_CpG)
  }
  
  if (!is.null(specific_pos)) {
    plot_df <- plot_df %>% filter(chr==specific_pos)
  }
  
  # show sample name of outliers
  outlier_data <- plot_df %>%
    group_by(!!rlang::sym(plotting_column), phenotype) %>%
    mutate(outliers = list(boxplot.stats(meth_perc)$out)) %>%
    filter(map_lgl(meth_perc, ~.x %in% outliers[[1]])) %>%
    select(-outliers)
  
  base_plot <- ggplot(plot_df, aes(x=plotting_x, y=meth_perc, fill=visualization_factor))
  
  # boxplot
  if (show_boxplot == TRUE) {
    if (!is.null(target_sample)) {
      base_plot <- base_plot + geom_boxplot(data = plot_df %>% filter(sample != target_sample), outlier.shape = NA)
    } else {
      base_plot <- base_plot + geom_boxplot(outlier.shape = NA)
    }
  }
  
  # points
  if (show_points == TRUE) {
    base_plot <- base_plot + geom_point(aes(color=visualization_factor), 
                                        position = position_jitterdodge(jitter.width = jitter_width, seed = jitter_seed), 
                                        alpha = 0.25)
  }
  
  # when plotting the target sample separately, make the points slightly larger
  if (!is.null(target_sample)) {
    base_plot <- base_plot + geom_point(data = plot_df %>% filter(sample == target_sample), 
                                        aes(color=visualization_factor), 
                                        shape=23,
                                        size=1.5)
  }
  
  # annotate outliers with sample labels
  base_plot <- base_plot + geom_text(data = outlier_data, aes(label = sample, color = visualization_factor), 
                                     position = position_jitterdodge(jitter.width = jitter_width, seed = jitter_seed), 
                                     hjust = 1.5, size = 3)
  
  return(base_plot + theme_classic() +
           guides(fill=guide_legend(title="Sample"), color="none") +
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
           theme(plot.title = element_text(hjust = 0.5)) +
           xlab("CpG position") + ylab("Percent methylation (%)") + ggtitle(input_gene))
}

print(plot_by_gene(CpG_df, "H19-BP1X", min_cov=25, show_boxplot=TRUE))

for (gene in target_sites) {
  print(gene)
  png(filename=paste("methylation-patterns-per-site/", gene, ".png", sep=""), width=1280, height=720)
  print(plot_by_gene(CpG_df, gene, min_cov=25, show_boxplot=TRUE, plot_chr = TRUE))
  dev.off();
}
