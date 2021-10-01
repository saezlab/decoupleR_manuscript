# List of the methods to call
stats_list = list(c('udt','mdt','aucell','wmean','wsum','ulm','mlm','viper','gsva','ora','fgsea'))

# Random seed
seed = 42
nproc = 4

# List of options for each method
opts_list <- list(list(
  aucell = list(nproc=nproc, seed=seed),
  fgsea = list(times=100, nproc=nproc, seed=seed),
  gsva = list(verbose = FALSE, method = "gsva"),
  mdt = list(sparse=FALSE, center=FALSE, na.rm=FALSE, trees = 10, min_n = 20, nproc = nproc, seed=seed),
  mlm = list(sparse=FALSE, center=FALSE, na.rm=FALSE),
  ora = list(n_up=300, n_bottom=300, n_background=20000, with_ties = TRUE),
  udt = list(sparse=FALSE, center=FALSE, na.rm=FALSE, min_n = 20, seed=seed),
  ulm = list(sparse=FALSE, center=FALSE, na.rm=FALSE),
  viper = list(verbose = FALSE, minsize = 0, pleiotropy = T, eset.filter = F),
  wmean = list(times=100, sparse=TRUE, randomize_type = "rows", seed=seed),
  wsum = list(times=100, sparse=TRUE, randomize_type = "rows", seed=seed)
))
