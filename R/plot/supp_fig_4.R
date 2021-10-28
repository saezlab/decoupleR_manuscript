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
get_auc_df <- function(df, .type){
  .type <- enquo(.type)
  df %>%
    select(set_name, statistic, !!.type) %>%
    mutate(!!.type := map(!!.type, function(df){
      df %>%
        group_by(run) %>%
        summarize(raw_auc = unique(raw_auc)) %>%
        pull(raw_auc)
    })) %>%
    unnest(cols = c(!!.type)) %>%
    filter(!statistic %in% c('aucell','ora','norm_fgsea','fgsea','gsva'))
}

get_auc_boxplot <- function(df, .type, ylabel='AUROC'){
  .type <- enquo(.type)

  order <- df %>%
    group_by(set_name, statistic) %>%
    summarize(median=median(!!.type), .groups='drop') %>%
    arrange(desc(median)) %>%
    distinct(statistic) %>%
    pull(statistic)
  df$statistic <- factor(df$statistic, levels = rev(order))

  ggplot(df, aes(x=statistic, y=!!.type, color=set_name)) +
    geom_boxplot() +
    xlab('Methods') +
    ylab(ylabel) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

get_auc_scatter <- function(df){
  min_lim <- floor(min(c(df$roc, df$prc)) * 100)/100
  max_lim <- ceiling(max(c(df$roc, df$prc)) * 100)/100
  ggplot(df, aes(x=roc, y=prc, label=statistic, color=set_name)) +
    geom_point() +
    geom_text_repel(max.overlaps = Inf, max.time=5, max.iter=1000000) +
    theme(text = element_text(size=14)) +
    xlab('AUROC') +
    ylab('AUPRC') +
    xlim(min_lim,max_lim) +
    ylim(min_lim,max_lim) +
    theme_bw()
}

# Read
rna_result <- readRDS(file.path('data', 'prc', 'rna_result.rds'))
php_result <- readRDS(file.path('data', 'prc', 'php_result.rds'))

# Generate data-frames
php_roc_df <- get_auc_df(php_result, roc)
php_prc_df <- get_auc_df(php_result, prc)

php_auc_df <- php_roc_df %>%
  left_join(php_prc_df) %>%
  group_by(set_name, statistic) %>%
  summarise(roc = median(roc), prc = median(prc), .groups='drop')

# Generate plots
php_roc_boxp <- get_auc_boxplot(php_roc_df, roc, 'AUROC')
php_prc_boxp <- get_auc_boxplot(php_prc_df, prc, 'AUPRC') + theme(legend.position="none")
php_auc_scatt <- get_auc_scatter(php_auc_df) + theme(legend.position="none")

# Save
pdf(file = file.path(path_figs, 'supp_fig_4.pdf'),
    width = 10, # The width of the plot in inches
    height = 5) # The height of the plot in inches
((php_roc_boxp / php_prc_boxp) | php_auc_scatt) +
  plot_layout(guides = 'collect', widths = c(1, 1))  +
  plot_annotation(tag_levels = 'A')
dev.off()

# Test significance better performance
df <- php_roc_df %>%
  left_join(php_prc_df) %>%
  pivot_longer(cols=c(roc, prc)) %>%
  select(-name) %>%
  filter(statistic != 'norm_wsum')
unweighted <- df %>%
  filter(set_name == 'unweighted') %>%
  pull(value)
weighted <- df %>%
  filter(set_name == 'weighted') %>%
  pull(value)
p_value <- wilcox.test(weighted, unweighted, alternative = "g")$p.value
print(paste0('p-value weighted methods have greater AUC: ', p_value))





