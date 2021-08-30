library(decoupleRBench)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)

# Benchmark data
bexpr_loc <- file.path('data', "dorothea_bench_expr.rds")
# benchmark metadata
bmeta_loc <- file.path('data', "dorothea_bench_meta.rds")

# List of the methods to call
stats_list = list(c('aucell','mean','pscira','scira','viper','gsva','ora','fgsea'))

# List of options for each method
opts_list <- list(list(
  aucell = list(nCores = 4),
  mean = list(times=100, .mor = "mor"),
  pscira = list(times=100, .mor = "mor"),
  scira = list(.mor = "mor", fast = FALSE, center=FALSE),
  viper = list(verbose = FALSE,
               minsize = 0,
               .mor = "mor",
               .likelihood = "likelihood"),
  gsva = list(verbose = FALSE, method = "gsva"),
  ora = list(n_up=300, n_bottom=300, n_background=20000),
  fgsea = list(force_ties = T, options = list(nproc=4))
))

# Net path
net_loc <- file.path('data', 'dorothea.rds')

# Design
design <- tibble(
  set_name = 'dorothea', # name of the set resource
  bench_name = "dbd", # name of the benchmark data
  stats_list = stats_list,
  opts_list = opts_list,
  bexpr_loc = bexpr_loc, # benchmark data location
  bmeta_loc = bmeta_loc, # metadata location
  source_loc = net_loc, # set source location
  source_col = "tf", # source name of the gene set source
  target_col = "target", # target name of the set source
  filter_col = "confidence", # column by which we wish to filter
  filter_crit = list(c('A','B','C')) # criteria by which we wish to filter
)

# Run benchmark
benchmark <- run_benchmark(
  .design = design, # provide input tibble
  .minsize = 10, # filter gene sets with size < 10
  .form = TRUE, # format the benchmark results
  .perform = TRUE, # evaluate benchmarking performance
  .silent = FALSE, # silently run the pipeline
  .downsample_pr = TRUE, # downsample TNs for precision-recall curve
  .downsample_roc = TRUE, # downsample TNs for ROC
  .downsample_times = 100, # downsampling iterations (not used here)
  .url_bool = FALSE # whether to load from url
)

# Correlation
df <- benchmark@bench_res %>%
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

cor_heat <- pheatmap(corr_matrix, cluster_rows = T,
                     cluster_cols = T)
idxs <- cor_heat$tree_row$order
corr_matrix <- corr_matrix[idxs,idxs]
corr_matrix[lower.tri(corr_matrix, diag=T)] <- NA


ComplexHeatmap::draw(lgd, x = unit(0.9, "npc"), y = unit(0.5, "npc"))

library(circlize)


cor_heat <- pheatmap(corr_matrix, color = colorRampPalette((RColorBrewer::brewer.pal(n = 7, name ="Greens")))(100),
                     display_numbers=F, cluster_rows = F, number_color='black', border_color=NA,
                     cluster_cols = F, na_col=NA, cellwidth = 15, cellheight = 15,
                     legend_breaks = c(0, 0.25, 0.50, 0.75, 1.0), legend=F,
                     show_rownames = T, show_colnames = T)



# Jaccard
stat_acts <- benchmark@bench_res %>%
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

jacc_matrix <- matrix(0, ncol(df), ncol(df))
colnames(jacc_matrix) <- colnames(df)
rownames(jacc_matrix) <- colnames(df)

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
    jacc_matrix[name_a,name_b] <- mean(jacs)
  }
}

jacc_matrix <- jacc_matrix[idxs,idxs]
jacc_matrix[upper.tri(jacc_matrix, diag=T)] <- NA
jac_heat <- pheatmap(jacc_matrix, color = colorRampPalette((RColorBrewer::brewer.pal(n = 7, name ="Reds")))(100),
                     display_numbers=F, cluster_rows = F, number_color='black', border_color=NA,
                     cluster_cols = F, na_col=NA, cellwidth = 15, cellheight = 15,
                     legend_breaks = c(0, 0.25, 0.50, 0.75, 1.0), legend = F,
                     show_rownames = T, show_colnames = T)




plot.new()
cor_heat
jac_heat
col_fun = colorRamp2(c(0, 0.25, 0.50, 0.75, 1.0),
                     colorRampPalette((RColorBrewer::brewer.pal(n = 7, name ="Greens")))(5))
lgd = ComplexHeatmap::Legend(col_fun = col_fun, title = "Corr")
ComplexHeatmap::draw(lgd, x = unit(0.19, "npc"), y = unit(0.6, "npc"))

col_fun = colorRamp2(c(0, 0.25, 0.50, 0.75, 1.0),
                     colorRampPalette((RColorBrewer::brewer.pal(n = 7, name ="Reds")))(5))
lgd = ComplexHeatmap::Legend(col_fun = col_fun, title = "Jacc")
ComplexHeatmap::draw(lgd, x = unit(0.1, "npc"), y = unit(0.6, "npc"))

# AUROC
aucs <- benchmark@bench_res %>%
  select(statistic, roc) %>%
  mutate(roc = map(roc, function(df){
    df %>%
      group_by(run) %>%
      summarize(raw_auc = unique(raw_auc)) %>%
      pull(raw_auc)
  })) %>%
  unnest(cols = c(roc))

ggplot(aucs, aes(x=forcats::fct_reorder(statistic, roc, .fun = median, .desc =TRUE), y=roc)) +
  theme_classic() +
  geom_boxplot() +
  ylim(0.5,1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# AUPR
prcs <- benchmark@bench_res %>%
  select(statistic, prc) %>%
  mutate(prc = map(prc, function(df){
    df %>%
      group_by(run) %>%
      summarize(raw_auc = unique(raw_auc)) %>%
      pull(raw_auc)
  })) %>%
  unnest(cols = c(prc))

ggplot(prcs, aes(x=forcats::fct_reorder(statistic, prc, .fun = median, .desc =TRUE), y=prc)) +
  theme_classic() +
  geom_boxplot() +
  ylim(0.5,1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Both
both <- aucs %>%
  left_join(prcs) %>%
  group_by(statistic) %>%
  summarise(roc = median(roc), prc = median(prc))

library(ggrepel)

ggplot(both, aes(x=roc, y=prc, label=statistic, color=statistic)) +
  theme_classic() +
  geom_point() +
  geom_text_repel(size=5)+
  theme(legend.position = "none") +
  ylim(0.5,0.7) +
  xlim(0.5,0.7)

