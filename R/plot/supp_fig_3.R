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

get_corr_df <- function(df){
  # Split df by normal, add and del
  lst_dfs <- df %>%
    select(mode, perm, perc, activity) %>%
    unnest(activity) %>%
    select(mode, perm, perc, statistic, source, id, score) %>%
    group_by(mode) %>%
    group_split()
  add_df <- lst_dfs[[1]]
  del_df <- lst_dfs[[2]]
  normal_df <- lst_dfs[[3]] %>%
    arrange(source, id, statistic) %>%
    pivot_wider(names_from=statistic, values_from=score) %>%
    unite(pair, source, id, sep='.')

  # Compute corrs
  map(list(add_df, del_df), function(noise_df){
    noise_df %>%
      group_by(perc) %>%
      group_split() %>%
      map(function(df){
        df %>%
          group_by(perm) %>%
          group_split() %>%
          map(function(df2){
            perc_n <- unique(df2$perc)
            mode_n <- unique(df2$mode)
            perm_n <- unique(df2$perm)
            tmp <- df2 %>%
              arrange(source, id, statistic) %>%
              pivot_wider(names_from=statistic, values_from=score) %>%
              unite(pair, source, id, sep='.') %>%
              filter(pair %in% intersect(.$pair, normal_df$pair)) %>%
              select(-mode, -perm, -perc)
            normal_df_tmp <- normal_df %>%
              filter(pair %in% intersect(.$pair, tmp$pair))
            tmp <- tmp %>%
              select(-pair)
            map(colnames(tmp), function(name){
              tibble(
                mode = mode_n,
                perc = perc_n,
                perm = perm_n,
                statistic = name,
                corr = cor(normal_df_tmp[[name]], tmp[[name]])
              )
            }) %>%
              bind_rows()
          }) %>%
          bind_rows()
      }) %>%
      bind_rows()
  }) %>%
    bind_rows()
}

get_corr_plot <- function(df, mode_noise, title){
  df %>%
    filter(mode == mode_noise) %>%
    ggplot(
      aes(x=forcats::fct_reorder(statistic, corr, .fun = median, .desc =TRUE),
          y=corr,
          color=perc)
    ) +
    theme_light() +
    geom_boxplot() +
    theme(text = element_text(size=14),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    xlab('Methods') +
    ylab('Correlation') +
    ggtitle(title) +
    ylim(0.5,1)
}

test_sign <- function(df){
  df <- df %>%
    group_by(mode) %>%
    group_split() %>%
    map(function(df){
      pull(df, corr)
    })
  wilcox.test(df[[1]], df[[2]], alternative = "g")$p.value
}

# Read
rna_result <- read_rds(file.path('data', 'prc', 'rna_noise.rds'))
php_result <- read_rds(file.path('data', 'prc', 'php_noise.rds'))

# Generate data-frames
rna_corr_df <- get_corr_df(rna_result)
php_corr_df <- get_corr_df(php_result)

# Test significance
median_corrs <- full_join(rna_corr_df, php_corr_df) %>%
  group_by(mode) %>%
  group_split() %>%
  map(function(df){
    df %>%
      group_by(mode) %>%
      summarize(median_corr = median(corr))
  }) %>%
  bind_rows()

pval <- test_sign(full_join(rna_corr_df, php_corr_df))
print(median_corrs)
print(paste0('wilcoxon pvalue: ', pval))

# Generate plots
rna_add_box <- get_corr_plot(rna_corr_df, mode='add', title='Addition')
rna_del_box <- get_corr_plot(rna_corr_df, mode='del', title='Deletion')
php_add_box <- get_corr_plot(php_corr_df, mode='add', title='Addition')
php_del_box <- get_corr_plot(php_corr_df, mode='del', title='Deletion')

# Merge together and save
pdf(file = file.path(path_figs, 'supp_fig_3.pdf'),
    width = 9, # The width of the plot in inches
    height = 9) # The height of the plot in inches
rna_add_box + rna_del_box + php_add_box + php_del_box +
  plot_layout(guides = 'collect')  +
  plot_annotation(tag_levels = 'A')
dev.off()



