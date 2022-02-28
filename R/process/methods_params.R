# List of the methods to call
stats_list = list(c('aucell','wmean','wsum','ulm','mlm','viper','gsva','ora','fgsea','udt','mdt'))

# Random seed
seed = 42
nproc = 8

# List of options for each method
opts_list <- list(list(
  udt = list(sparse=FALSE, center=FALSE, na.rm=FALSE, min_n = 20, seed=seed),
  mdt = list(sparse=FALSE, center=FALSE, na.rm=FALSE, trees = 10, min_n = 20,
             nproc = nproc, seed=seed),
  aucell = list(nproc=nproc, seed=seed),
  wmean = list(times=100, sparse=FALSE, randomize_type = "rows", seed=seed),
  ulm = list(sparse=FALSE, center=FALSE, na.rm=FALSE),
  mlm = list(sparse=FALSE, center=FALSE, na.rm=FALSE),
  viper = list(verbose=FALSE, pleiotropy=T, eset.filter=F),
  gsva = list(verbose = FALSE, method = "gsva"),
  ora = list(n_up=300, n_bottom=300, n_background=20000, with_ties = TRUE, seed=seed),
  fgsea = list(times=100, nproc=nproc, seed=seed),
  wsum = list(times=100, sparse=FALSE, randomize_type = "rows", seed=seed)
))
