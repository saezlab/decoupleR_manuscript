library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(ggrepel)
library(patchwork)
library(ggplotify)

# Create dir
path_figs <- file.path('figures')
dir.create(path_figs, showWarnings = F, recursive = T)

# Plot functions
get_df_scores <- function(df){
  df %>%
    select(activity) %>%
    unnest(activity) %>%
    select(statistic, score)
}

get_dist_plot <- function(df){
  df %>%
    ggplot(aes(x=score, color=statistic, fill=statistic)) +
    geom_density() +
    theme(text = element_text(size=12)) +
    facet_wrap(~ statistic, scales='free') +
    xlab('scores') +
    theme_bw()
}

# Read
rna_result <- readRDS(file.path('data', 'prc', 'rna_result.rds'))
php_result <- readRDS(file.path('data', 'prc', 'php_result.rds'))

# Generate dfs
rna_scores <- get_df_scores(rna_result) %>% mutate(dataset = 'transcriptomic')
php_scores <- get_df_scores(php_result) %>% mutate(dataset = 'phosphoproteomic')

# Generate plots
get_dist_plot <- function(df){
  df %>%
    ggplot(aes(x=score, color=statistic, fill=statistic)) +
    geom_density() +
    theme(text = element_text(size=12)) +
    facet_wrap(~ statistic, scales='free') +
    xlab('scores') +
    ylab('densities') +
    theme_bw() +
    theme(legend.position = "none")
}

# Merge together and save
pdf(file = file.path(path_figs, 'supp_fig_1.pdf'),
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
(get_dist_plot(rna_scores)) / get_dist_plot(php_scores) +
  plot_layout(guides = 'collect')  +
  plot_annotation(tag_levels = 'A')
dev.off()


