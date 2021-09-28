# List of the methods to call
stats_list = list(c('udt','mdt','aucell','wmean','wsum','ulm','mlm','viper','gsva','ora','fgsea'))

# List of options for each method
opts_list <- list(list(
  udt = list(min_n = 20),
  mdt = list(trees = 1000, min_n = 20, nproc = 4),
  aucell = list(nproc=4),
  wmean = list(times=1000),
  wsum = list(times=1000),
  ulm = list(center=FALSE),
  mlm = list(center=FALSE),
  viper = list(verbose = FALSE,
               minsize = 0,
               .mor = "mor",
               .likelihood = "likelihood",
               pleiotropy = T,
               eset.filter = F),
  gsva = list(verbose = FALSE, method = "gsva"),
  ora = list(n_up=300, n_bottom=300, n_background=20000),
  fgsea = list(force_ties = T, times=1000, nproc=4)
))
