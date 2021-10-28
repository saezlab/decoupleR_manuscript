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
get_time_df <- function(df){
  df %>%
    mutate(Time = map(activity, function(df){
      df$statistic_time[[1]]
    })) %>%
    mutate(Time = as.double(unlist(Time))) %>%
    select(statistic, Time)
}

get_time_plot <- function(df){
  ggplot(df, aes(x=forcats::fct_reorder(statistic, Time, .fun = median, .desc =TRUE),y=Time)) +
    geom_boxplot(outlier.size=0) +
    xlab('Methods') +
    ylab('Time (s)') +
    scale_y_continuous(
      breaks = seq(0,
                   ceiling(max(df$Time[!is.na(df$Time)]) / 30) * 30, by = 30)
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

# Read
rna_result <- readRDS(file.path('data', 'prc', 'rna_noise.rds'))
php_result <- readRDS(file.path('data', 'prc', 'php_noise.rds'))

# Generate data-frames
rna_time_df <- get_time_df(rna_result)
php_time_df <- get_time_df(php_result)

# Generate figures
rna_time_p <- get_time_plot(rna_time_df)
php_time_p <- get_time_plot(php_time_df)

# Merge together and save
pdf(file = file.path(path_figs, 'supp_fig_7.pdf'),
    width = 6, # The width of the plot in inches
    height = 3) # The height of the plot in inches
rna_time_p + php_time_p +
  plot_annotation(tag_levels = 'A')
dev.off()

# Report median speed
top_performers <- read.csv(file.path(path_figs, 'supp_tab_2.csv')) %>%
  filter(p_value < 0.05) %>%
  pull(statistic)

get_sample_n <- function(df){
  df[1,] %>%
    select(activity) %>%
    unnest(activity) %>%
    pull(id) %>%
    unique() %>%
    length()
}

get_regult_n <- function(df){
  df[1,] %>%
    select(activity) %>%
    unnest(activity) %>%
    pull(source) %>%
    unique() %>%
    length()
}

rna_sample_n <- get_sample_n(rna_result)
rna_regult_n <- get_regult_n(rna_result)
php_sample_n <- get_sample_n(php_result)
php_regult_n <- get_regult_n(php_result)


rna_speed_df <- rna_time_df %>%
  filter(statistic %in% top_performers) %>%
  mutate(Speed = Time / rna_sample_n / rna_regult_n) %>%
  select(statistic, Speed)

php_speed_df <- php_time_df %>%
  filter(statistic %in% top_performers) %>%
  mutate(Speed = Time / php_sample_n / php_regult_n) %>%
  select(statistic, Speed)

median_speed <- rna_speed_df %>%
  bind_rows(php_speed_df) %>%
  summarise(median=median(Speed)) %>%
  pull(median)

print(paste0('Median speed of top performers: ',
             formatC(median_speed, format = "e", digits = 2)))

