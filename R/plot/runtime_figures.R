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

# RNA
rna_noise <- readRDS("~/saezlab/decoupleR_manuscript/data/prc/rna_noise.rds")

rna_time <- rna_noise %>%
  mutate(Time = map(activity, function(df){
    df$statistic_time[[1]]
  })) %>%
  mutate(Time = as.double(unlist(Time))) %>%
  select(statistic, Time)

p_rna <- ggplot(rna_time,
       aes(x=forcats::fct_reorder(statistic, Time, .fun = median, .desc =TRUE),
           y=Time)
) +
  theme_light() +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab('Methods') +
  ylab('Time (s)') +
  scale_y_continuous(
    breaks = seq(0,
                 ceiling(max(rna_time$Time[!is.na(rna_time$Time)]) / 30) * 30, by = 30)
    ) +
  ggtitle('Transcriptomics data-set')


# PHP
php_noise <- readRDS("~/saezlab/decoupleR_manuscript/data/prc/php_noise.rds")

php_time <- php_noise %>%
  mutate(Time = map(activity, function(df){
    df$statistic_time[[1]]
  })) %>%
  mutate(Time = as.double(unlist(Time))) %>%
  select(statistic, Time)

p_php <- ggplot(php_time,
       aes(x=forcats::fct_reorder(statistic, Time, .fun = median, .desc =TRUE),
           y=Time)
) +
  theme_light() +
  geom_boxplot() +
  ggtitle('Phospho-protemoics data-set') +
  xlab('Methods') +
  ylab('Time (s)') +
  scale_y_continuous(
    breaks = seq(0,
                 ceiling(max(php_time$Time[!is.na(php_time$Time)]) / 15) * 15, by = 15)
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Write
pdf(file = file.path(path_figs, 'runtime.pdf'),
    width = (1*4), # The width of the plot in inches
    height = (2*4)) # The height of the plot in inches
p_rna / p_php + plot_annotation(tag_levels = 'a') +
  theme(plot.tag = element_text(face = 'bold'))
dev.off()
