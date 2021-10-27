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

get_auc_df <- function(df, .type){
  .type <- enquo(.type)
  df %>%
    select(mode, perm, perc, statistic, !!.type) %>%
    unnest(!!.type) %>%
    mutate(mode = factor(mode, levels=c('normal', 'add', 'del')))
}

get_auc_boxplot <- function(df){
  ggplot(df, aes(x=mode, y=raw_auc, color=perc)) +
    geom_boxplot() +
    theme(text = element_text(20),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    xlab('') +
    ylab('AUCs') +
    facet_wrap(~ statistic) +
    theme_bw()
}

# Read
rna_result <- read_rds(file.path('data', 'prc', 'rna_noise.rds'))
php_result <- read_rds(file.path('data', 'prc', 'php_noise.rds'))

# Process dfs
rna_roc <- get_auc_df(rna_result, roc)
rna_prc <- get_auc_df(rna_result, prc)
php_roc <- get_auc_df(php_result, roc)
php_prc <- get_auc_df(php_result, prc)
df <- full_join(full_join(rna_roc, rna_prc), full_join(php_roc, php_prc))

# Plot
plt <- get_auc_boxplot(df)

# Test sign
two.way <- aov(raw_auc ~ mode, data = mutate(df,mode = factor(
  mode, levels=c('del', 'add', 'normal'))))

summary(two.way)
tukey.two.way<-TukeyHSD(two.way)

tukey.two.way

# Save
pdf(file = file.path(path_figs, 'supp_fig_6.pdf'),
    width = 7, # The width of the plot in inches
    height = 6) # The height of the plot in inches
plt
dev.off()

