library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(ggrepel)
library(patchwork)

# Create dir
path_figs <- file.path('figures')
dir.create(path_figs, showWarnings = F, recursive = T)

# Read
php_fname <- file.path('data', 'prc', 'php_noise.rds')
php_result <- readRDS(php_fname) %>%
  tidyr::separate(set_name, c('mode','perm'), sep='(?=[[:digit:]]+)') %>%
  rename('perc'=bench_name) %>%
  mutate(perc = paste0(as.double(perc)*100, ' %'))

###
# Correlations of activities by noise
###
# Divide by normal, add and del
lst_dfs <- php_result %>%
  select(mode, perm, perc, activity) %>%
  unnest(activity) %>%
  select(mode, perm, perc, statistic, tf, id, score) %>%
  group_by(mode) %>%
  group_split()
add_df <- lst_dfs[[1]]
del_df <- lst_dfs[[2]]
normal_df <- lst_dfs[[3]] %>%
  arrange(tf, id, statistic) %>%
  pivot_wider(names_from=statistic, values_from=score) %>%
  unite(pair, tf, id, sep='.')

# Correlations
corrs <- map(list(add_df, del_df), function(noise_df){
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
            arrange(tf, id, statistic) %>%
            pivot_wider(names_from=statistic, values_from=score) %>%
            unite(pair, tf, id, sep='.') %>%
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

# Plots
n_max <- max(corrs$corr)
n_min <- min(corrs$corr)
n_max <- n_max + n_max * 0.05
n_min <- n_min + n_min * 0.05
add_cors_p <- corrs %>%
  filter(mode=='add') %>%
  ggplot(
    aes(x=forcats::fct_reorder(statistic, corr, .fun = median, .desc =TRUE),
        y=corr,
        color=perc)
  ) +
  theme_classic() +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab('Methods') +
  ylab('Correlation') +
  theme(aspect.ratio=1) +
  ggtitle('Adding random edges') +
  ylim(n_min, n_max)
del_cors_p <- corrs %>%
  filter(mode=='del') %>%
  ggplot(
    aes(x=forcats::fct_reorder(statistic, corr, .fun = median, .desc =TRUE),
        y=corr,
        color=perc)
  ) +
  theme_classic() +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab('Methods') +
  ylab('Correlation') +
  theme(aspect.ratio=1) +
  ggtitle('Deleting random edges') +
  ylim(n_min, n_max)

# Write
pdf(file = file.path(path_figs, 'php_noise_corr.pdf'),
    width = 14, # The width of the plot in inches
    height = 7) # The height of the plot in inches
add_cors_p + del_cors_p
dev.off()

###
# Difference of AUROC/AUPRC by noise
###
# AUROC
aucs <- php_result %>%
  select(-prc, -activity) %>%
  mutate(roc = map(roc, function(df){
    df %>%
      group_by(run) %>%
      summarize(raw_auc = unique(raw_auc)) %>%
      pull(raw_auc)
  })) %>%
  unnest(cols = c(roc))

# AUPRC
prcs <- php_result %>%
  select(-roc, -activity) %>%
  mutate(prc = map(prc, function(df){
    df %>%
      group_by(run) %>%
      summarize(raw_auc = unique(raw_auc)) %>%
      pull(raw_auc)
  })) %>%
  unnest(cols = c(prc))

# Both
both <- aucs %>%
  left_join(prcs) %>%
  group_by(mode, perc, statistic) %>%
  summarise(roc = median(roc), prc = median(prc), .groups='drop')

# Differences
normal <- both %>%
  filter(mode == 'normal')
noise <- both %>%
  filter(mode != 'normal') %>%
  group_by(statistic) %>%
  group_split() %>%
  map(function(df){
    stat_name <- df$statistic[[1]]
    og_df <- normal %>%
      filter(statistic == stat_name)
    og_roc <- og_df %>% pull(roc)
    og_prc <- og_df %>% pull(prc)
    df %>%
      mutate(
        roc = roc - og_roc,
        prc = prc - og_prc
      )
  }) %>%
  bind_rows()

# Plots
diff_add_p <- ggplot(dplyr::filter(noise, mode=='add'),
       aes(x=roc, y=prc, color=statistic, size=perc)) +
  theme_classic() +
  geom_point() +
  xlab('ΔAUROC') +
  ylab('ΔAUPRC') +
  theme(aspect.ratio=1) +
  scale_color_brewer(palette = "Paired") +
  geom_hline(yintercept = 0, linetype = 'dashed', alpha=0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed', alpha=0.5) +
  ggtitle("Adding random edges")
diff_del_p <- ggplot(dplyr::filter(noise, mode=='del'),
       aes(x=roc, y=prc, color=statistic, size=perc)) +
  theme_classic() +
  geom_point() +
  xlab('ΔAUROC') +
  ylab('ΔAUPRC') +
  theme(aspect.ratio=1) +
  scale_color_brewer(palette = "Paired") +
  geom_hline(yintercept = 0, linetype = 'dashed', alpha=0.5) +
  geom_vline(xintercept = 0, linetype = 'dashed', alpha=0.5) +
  ggtitle("Deleting random edges")
# Write
pdf(file = file.path(path_figs, 'php_noise_diff_auc.pdf'),
    width = 14, # The width of the plot in inches
    height = 7,
    encoding='Greek'
    ) # The height of the plot in inches
diff_add_p + diff_del_p
dev.off()



