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
  statistics <- unique(df$statistic)
  lst_dfs <- df %>%
    select(mode, perm, perc, activity) %>%
    unnest(activity) %>%
    select(mode, perm, perc, statistic, source, id, score) %>%
    mutate(pair = paste0(source, '.', id)) %>%
    group_by(mode) %>%
    group_split()
  unq_pairs <- lst_dfs[[2]] %>%
    group_by(perm, perc) %>%
    group_split() %>%
    map(function(df){
      df %>%
        pull(pair)
    }) %>%
    Reduce(intersect, .)
  add_df <- lst_dfs[[1]] %>%
    filter(pair %in% unq_pairs)
  del_df <- lst_dfs[[2]] %>%
    filter(pair %in% unq_pairs)
  nrm_df <- lst_dfs[[3]] %>%
    filter(pair %in% unq_pairs) %>%
    pivot_wider(names_from=statistic, values_from=score)

  map(list(add_df, del_df), function(noise_df){
    noise_df %>%
      group_by(perc, perm) %>%
      group_split() %>%
      map(function(df){
        perc_n <- unique(df$perc)
        mode_n <- unique(df$mode)
        perm_n <- unique(df$perm)
        df <- df %>%
          pivot_wider(names_from=statistic, values_from=score)
        map(statistics, function(statistic){
          tibble(
            mode = mode_n,
            perc = perc_n,
            perm = perm_n,
            statistic = statistic,
            corr = cor(df[[statistic]], nrm_df[[statistic]], method='spearman')
          )
        }) %>%
          bind_rows()
      }) %>%
      bind_rows()
    }) %>%
    bind_rows()
}

get_corr_plot <- function(df, mode_noise, title, min_corr=0){
  df %>%
    filter(mode == mode_noise) %>%
    ggplot(
      aes(x=forcats::fct_reorder(statistic, corr, .fun = median, .desc =TRUE),
          y=corr,
          color=perc)
    ) +
    geom_boxplot(outlier.size=0) +
    xlab('Methods') +
    ylab('Correlation') +
    ggtitle(title) +
    ylim(min_corr,1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

# Read
rna_result <- read_rds(file.path('data', 'prc', 'rna_noise.rds'))
php_result <- read_rds(file.path('data', 'prc', 'php_noise.rds'))

# Generate data-frames
rna_corr_df <- get_corr_df(rna_result)
php_corr_df <- get_corr_df(php_result)

# Generate plots
both_corr_df <- rbind(rna_corr_df, php_corr_df)
min_corr <- both_corr_df %>% pull(corr) %>% min()
min_corr <- min_corr - (min_corr * 0.05)
rna_add_box <- get_corr_plot(rna_corr_df, mode='add', title='Addition', min_corr=min_corr)
rna_del_box <- get_corr_plot(rna_corr_df, mode='del', title='Deletion', min_corr=min_corr)
php_add_box <- get_corr_plot(php_corr_df, mode='add', title='Addition', min_corr=min_corr)
php_del_box <- get_corr_plot(php_corr_df, mode='del', title='Deletion', min_corr=min_corr)

# Test sign
median_corr <- rna_corr_df %>%
  bind_rows(php_corr_df) %>%
  group_by(mode) %>%
  summarize(median_corr=median(corr))
print(median_corr)
add <- filter(both_corr_df, mode=='add')$corr
del <- filter(both_corr_df, mode=='del')$corr
test <- wilcox.test(add, del, alternative = "g", paired=T)
p_value <- formatC(test$p.value, format = "e", digits = 2)
W <- formatC(unname(test$statistic), format = "e", digits = 2)
N <- formatC(length(add), format = "e", digits = 2)
print(paste0('add-del: p_value=', p_value, '; W=', W,'; N=', N))

# Merge together and save
pdf(file = file.path(path_figs, 'supp_fig_5.pdf'),
    width = 6, # The width of the plot in inches
    height = 6) # The height of the plot in inches
rna_add_box + rna_del_box + php_add_box + php_del_box +
  plot_layout(guides = 'collect')  +
  plot_annotation(tag_levels = 'A')
dev.off()
