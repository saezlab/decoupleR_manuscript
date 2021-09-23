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
read_rds <- function(path){
  readRDS(path) %>%
    tidyr::separate(set_name, c('mode','perm'), sep='(?=[[:digit:]]+)') %>%
    rename('perc'=bench_name) %>%
    mutate(perc = paste0(as.double(perc)*100, ' %'))
}

# Read
rna_result <- read_rds(file.path('data', 'prc', 'rna_noise.rds'))
php_result <- read_rds(file.path('data', 'prc', 'php_noise.rds'))

# Plot functions
get_auc_df <- function(df, .type){
  .type <- enquo(.type)
  df %>%
    select(mode, perm, perc, statistic, !!.type) %>%
    unnest(!!.type)
}

get_auc_df(rna_result, roc)

rocs <- rna_result %>%
  select(-prc, -activity) %>%
  mutate(roc = map(roc, function(df){
    df %>%
      group_by(run) %>%
      summarize(raw_auc = unique(raw_auc)) %>%
      pull(raw_auc)
  })) %>%
  unnest(cols = c(roc))

# AUPRC
prcs <- rna_result %>%
  select(-roc, -activity) %>%
  mutate(prc = map(prc, function(df){
    df %>%
      group_by(run) %>%
      summarize(raw_auc = unique(raw_auc)) %>%
      pull(raw_auc)
  })) %>%
  unnest(cols = c(prc))

df <- rna_result %>%
  full_join(php_result) %>%
  select(-prc, -activity) %>%
  mutate(roc = map(roc, function(df){
    df %>%
      group_by(run) %>%
      summarize(raw_auc = unique(raw_auc)) %>%
      pull(raw_auc)
  })) %>%
  unnest(cols = c(roc))

two.way <- aov(roc ~ mode, data = df)

summary(two.way)
tukey.two.way<-TukeyHSD(two.way)

tukey.two.way



