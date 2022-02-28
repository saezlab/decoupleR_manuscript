library(decoupleR)
library(tidyverse)

# Get prc path
raw_path <- file.path('data', 'raw')
prc_path <- file.path('data', 'prc')
dir.create(prc_path, showWarnings = F, recursive = T)

# Get toy data-set
lst <- decoupleR::get_toy_data()
mat <- lst$mat
network <- lst$network

write.csv(mat, file.path(raw_path, 'toy_mat.csv'))
write.csv(network, file.path(raw_path, 'toy_net.csv'), row.names = F)

# Args
seed = 42
nproc = 8
opts_list <- list(
  udt = list(sparse=FALSE, center=FALSE, na.rm=FALSE, min_n = 3, seed=seed),
  mdt = list(sparse=FALSE, center=FALSE, na.rm=FALSE, trees = 10, min_n = 3,
             nproc = nproc, seed=seed),
  aucell = list(nproc=nproc, seed=seed, aucMaxRank=3),
  wmean = list(times=100, sparse=FALSE, randomize_type = "rows", seed=seed),
  wsum = list(times=100, sparse=FALSE, randomize_type = "rows", seed=seed),
  ulm = list(sparse=FALSE, center=FALSE, na.rm=FALSE),
  mlm = list(sparse=FALSE, center=FALSE, na.rm=FALSE),
  viper = list(verbose=FALSE, pleiotropy=T, eset.filter=F),
  gsva = list(verbose = FALSE, method = "gsva"),
  ora = list(n_up=3, n_bottom=0, n_background=20000, with_ties = TRUE, seed=seed),
  fgsea = list(times=100, nproc=nproc, seed=seed)
)

# Run
res <- decoupleR::decouple(mat,
                           network,
                           args = opts_list,
                           statistics = 'all',
                           minsize=0) %>%
  select(-run_id)

# Rename
res$statistic[res$statistic=='corr_wmean'] <- 'wmean_corr'
res$statistic[res$statistic=='norm_wmean'] <- 'wmean_norm'
res$statistic[res$statistic=='corr_wsum'] <- 'wsum_corr'
res$statistic[res$statistic=='norm_wsum'] <- 'wsum_norm'
res$statistic[res$statistic=='fgsea'] <- 'gsea'
res$statistic[res$statistic=='norm_fgsea'] <- 'gsea_norm'
colnames(res) <- c("method", "source", "sample", "score", "pval")

# Sort, filter
res <- res %>%
  arrange(method, source, sample) %>%
  select(sample, source, score, method, pval) %>%
  mutate(score=ifelse(is.finite(score), score, NA))

# Store
write.csv(res, file.path(prc_path, 'toy_acts_R.csv'), row.names = F)
