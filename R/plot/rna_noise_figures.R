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
rna_fname <- file.path('data', 'prc', 'rna_noise.rds')
rna_result <- readRDS(rna_fname) %>%
  tidyr::separate(set_name, c('mode','perm'), sep='(?=[[:digit:]]+)') %>%
  rename('perc'=bench_name) %>%
  mutate(perc = paste0(as.double(perc)*100, ' %'))

###
# Correlations of activities by noise
###
# Divide by normal, add and del
lst_dfs <- rna_result %>%
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
  theme_light() +
  geom_boxplot() +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
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
  theme_light() +
  geom_boxplot() +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab('Methods') +
  ylab('Correlation') +
  theme(aspect.ratio=1) +
  ggtitle('Deleting random edges') +
  ylim(n_min, n_max)

# Write
pdf(file = file.path(path_figs, 'rna_noise_corr.pdf'),
    width = 14, # The width of the plot in inches
    height = 7) # The height of the plot in inches
add_cors_p + del_cors_p + plot_annotation(tag_levels = 'a')
dev.off()

###
# Difference of AUROC/AUPRC by noise
###
# AUROC
aucs <- rna_result %>%
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

# Facet plot aucs
medians <- aucs %>%
  filter(mode == 'normal') %>%
  group_by(statistic) %>%
  summarise(med = median(roc)) %>%
  arrange(desc(med))
medians$statistic <- factor(medians$statistic, levels = medians$statistic)

p_roc_add <- aucs %>%
  filter(mode == 'add') %>%
  mutate(statistic=factor(statistic, levels=pull(arrange(medians, desc(med)),statistic))) %>%
  arrange(statistic) %>%
  ggplot(aes(x=perc, y=roc)) +
  geom_boxplot() +
  geom_hline(data = medians, aes(yintercept = med), linetype="dashed") +
  facet_grid(~statistic) +
  theme_light() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  ylim(0.5, 1) +
  xlab('') +
  ylab('AUROC') +
  ggtitle('Adding random edges')

p_roc_del <- aucs %>%
  filter(mode == 'del') %>%
  mutate(statistic=factor(statistic, levels=pull(arrange(medians, desc(med)),statistic))) %>%
  arrange(statistic) %>%
  ggplot(aes(x=perc, y=roc)) +
  geom_boxplot() +
  geom_hline(data = medians, aes(yintercept = med), linetype="dashed") +
  facet_grid(~statistic) +
  theme_light() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  ylim(0.5, 1) +
  xlab('') +
  ylab('AUROC') +
  ggtitle('Deleting random edges')

# Facet plots prcs
medians <- prcs %>%
  filter(mode == 'normal') %>%
  group_by(statistic) %>%
  summarise(med = median(prc)) %>%
  arrange(desc(med))
medians$statistic <- factor(medians$statistic, levels = medians$statistic)

p_prc_add <- prcs %>%
  filter(mode == 'add') %>%
  mutate(statistic=factor(statistic, levels=pull(arrange(medians, desc(med)),statistic))) %>%
  arrange(statistic) %>%
  ggplot(aes(x=perc, y=prc)) +
  geom_boxplot() +
  geom_hline(data = medians, aes(yintercept = med), linetype="dashed") +
  facet_grid(~statistic) +
  theme_light() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  ylim(0.5, 1) +
  xlab('') +
  ylab('AUPRC') +
  ggtitle('Adding random edges')

p_prc_del <- prcs %>%
  filter(mode == 'del') %>%
  mutate(statistic=factor(statistic, levels=pull(arrange(medians, desc(med)),statistic))) %>%
  arrange(statistic) %>%
  ggplot(aes(x=perc, y=prc)) +
  geom_boxplot() +
  geom_hline(data = medians, aes(yintercept = med), linetype="dashed") +
  facet_grid(~statistic) +
  theme_light() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  ) +
  ylim(0.5, 1) +
  xlab('') +
  ylab('AUPRC') +
  ggtitle('Deleting random edges')

pdf(file = file.path(path_figs, 'rna_noise_auc.pdf'),
    width = (4*4), # The width of the plot in inches
    height = (4*4)) # The height of the plot in inches
p <- p_roc_add / p_roc_del / p_prc_add / p_prc_del
p + plot_annotation(tag_levels = 'a')
dev.off()

