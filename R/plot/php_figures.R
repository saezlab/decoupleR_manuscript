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

# Read
php_fname <- file.path('data', 'prc', 'php_result.rds')
php_result <- readRDS(php_fname)

###
# Similarties between activities
###
# Correlation
df <- php_result %>%
  dplyr::filter(set_name == 'weighted') %>%
  dplyr::select(set_name, filter_crit, activity) %>%
  tidyr::unnest(activity) %>%
  select(set_name, filter_crit, id, tf, statistic, score) %>%
  pivot_wider(names_from=statistic, values_from=score) %>%
  select(-set_name, -filter_crit, -id, -tf)
corr_matrix <- matrix(0, ncol(df), ncol(df))
colnames(corr_matrix) <- colnames(df)
rownames(corr_matrix) <- colnames(df)
for (name_a in colnames(df)) {
  for (name_b in colnames(df)) {
    corr_matrix[name_a,name_b] <- cor(abs(df[[name_a]]), abs(df[[name_b]]))
  }
}

# Get order of clustered methods
cor_heat <- pheatmap(corr_matrix, cluster_rows = T,
                     cluster_cols = T)
idxs <- cor_heat$tree_row$order
corr_matrix <- corr_matrix[idxs,idxs]

# Plot
celldim <- 10
cor_heat <- pheatmap(corr_matrix, color = colorRampPalette((RColorBrewer::brewer.pal(n = 7, name ="Greens")))(100),
                     display_numbers=F, cluster_rows = F, number_color='black', border_color=NA,
                     cluster_cols = T, na_col=NA, cellwidth = celldim, cellheight = celldim,
                     legend_breaks = c(0, 0.25, 0.50, 0.75, 1.0), legend=T,
                     show_rownames = T, show_colnames = T, main='Pearson correlation',
                     width = 4, height=4
                     )

# Jaccard
stat_acts <- php_result %>%
  dplyr::filter(set_name == 'weighted') %>%
  select(activity) %>%
  unnest(cols=c(activity)) %>%
  select(statistic, tf, id, score)

n_top <- length(stat_acts$tf %>% unique())
n_top <- ceiling(n_top * 0.05)
stat_acts <- stat_acts %>%
  group_by(statistic, id) %>%
  arrange(desc(abs(score))) %>%
  slice_head(n=n_top) %>%
  group_by(statistic, id) %>%
  select(-score) %>%
  nest(data=c(tf)) %>%
  pivot_wider(id_cols = id, names_from = statistic, values_from = data) %>%
  column_to_rownames('id')

jacc_matrix <- matrix(0, ncol(stat_acts), ncol(stat_acts))
colnames(jacc_matrix) <- colnames(stat_acts)
rownames(jacc_matrix) <- colnames(stat_acts)

jacc_idx <- function(a,b){
  n_inter <- length(intersect(a,b))
  n_union <- length(union(a,b))
  n_inter / n_union
}

for (name_a in colnames(stat_acts)) {
  for (name_b in colnames(stat_acts)) {
    jacs <- map_dbl(rownames(stat_acts), function(sample){
      tfs_a <- stat_acts[sample,name_a][[1]]$tf
      tfs_b <- stat_acts[sample,name_b][[1]]$tf
      jacc_idx(tfs_a, tfs_b)
    })
    jacc_matrix[name_a,name_b] <- round(mean(jacs), 4)
  }
}

# Get order of clustered methods
jac_heat <- pheatmap(jacc_matrix, cluster_rows = T,
                     cluster_cols = T, display_numbers = T)
idxs <- jac_heat$tree_row$order
jacc_matrix <- jacc_matrix[idxs,idxs]
jac_heat <-
pheatmap(jacc_matrix, color = colorRampPalette((RColorBrewer::brewer.pal(n = 7, name ="Reds")))(100),
                     display_numbers=F, cluster_rows = F, number_color='black', border_color=NA,
                     cluster_cols = T, na_col=NA, cellwidth = celldim, cellheight = celldim,
                     legend_breaks = c(0, 0.25, 0.50, 0.75, 1.0), legend=T,
                     show_rownames = T, show_colnames = T, main='Jaccard index',
                     width = 4, height=4
                     )

# Relationship between corr and jaccard
stats <- rownames(corr_matrix)
corr_matrix <- corr_matrix[stats,stats]
jacc_matrix <- jacc_matrix[stats,stats]
corr_jacc <- tibble(
  corr = corr_matrix[upper.tri(corr_matrix)],
  jacc = jacc_matrix[upper.tri(jacc_matrix)]
)
corr_jacc_p <- ggplot(corr_jacc, aes(x=corr,
                      y=jacc)) +
  theme_classic() +
  geom_point() +
  xlab('Correlations') +
  ylab('Jaccard Indexes') + theme(aspect.ratio=1)

# Write
pdf(file = file.path(path_figs, 'php_corr_jacc.pdf'),
    width = 12, # The width of the plot in inches
    height = 4.5) # The height of the plot in inches
as.ggplot(cor_heat) + as.ggplot(jac_heat) + corr_jacc_p
dev.off()

###
# Performance
###
# AUROC
aucs <- php_result %>%
  select(set_name, statistic, roc) %>%
  mutate(roc = map(roc, function(df){
    df %>%
      group_by(run) %>%
      summarize(raw_auc = unique(raw_auc)) %>%
      pull(raw_auc)
  })) %>%
  unnest(cols = c(roc)) %>%
  rename('Type' = set_name)

# AUPR
prcs <- php_result %>%
  select(set_name, statistic, prc) %>%
  mutate(prc = map(prc, function(df){
    df %>%
      group_by(run) %>%
      summarize(raw_auc = unique(raw_auc)) %>%
      pull(raw_auc)
  })) %>%
  unnest(cols = c(prc)) %>%
  rename('Type' = set_name)

# Both
both <- aucs %>%
  left_join(prcs) %>%
  group_by(Type, statistic) %>%
  summarise(roc = median(roc), prc = median(prc), .groups='drop')

# Plots
roc_p <- ggplot(aucs,
       aes(x=forcats::fct_reorder(statistic, roc, .fun = median, .desc =TRUE),
           y=roc,
           color=Type)
       ) +
  theme_classic() +
  geom_boxplot() +
  ylim(0.5,1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab('Methods') +
  ylab('AUROC')

prc_p <- ggplot(prcs, aes(x=forcats::fct_reorder(statistic, prc, .fun = median, .desc =TRUE),
                          y=prc, color=Type)) +
  theme_classic() +
  geom_boxplot() +
  ylim(0.5,1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab('Methods') +
  ylab('AUPRC')

both_p <- ggplot(both, aes(x=roc, y=prc, label=statistic, color=Type)) +
  theme_classic() +
  geom_point() +
  geom_text_repel() +
  xlab('AUROC') +
  ylab('AUPRC') +
  geom_line(aes(group = statistic), color='gray',
            arrow = arrow(length=unit(0.30,"cm"), type = "closed"))

# Write
pdf(file = file.path(path_figs, 'php_roc_prc.pdf'),
    width = (2*4), # The width of the plot in inches
    height = (3*4)) # The height of the plot in inches
(roc_p + prc_p) / both_p
dev.off()
