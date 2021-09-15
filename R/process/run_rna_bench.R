library(decoupleRBench)
library(dplyr)
library(tibble)

# Define data, metadata and network path
raw_path <- file.path('data', 'raw')
prc_path <- file.path('data', 'prc')
dir.create(prc_path, showWarnings = F, recursive = T)
expr_fname <- file.path(raw_path, "rna_expr.rds")
meta_fname <- file.path(raw_path, "rna_meta.rds")
netw_fname <- file.path(raw_path, 'dorothea.rds')

# List of the methods to call
stats_list = list(c('aucell','wmean','wsum','scira','viper','gsva','ora','fgsea'))

# List of options for each method
opts_list <- list(list(
  aucell = list(nCores=4),
  wmean = list(times=100, .mor = "mor"),
  wsum = list(times=100, .mor = "mor"),
  scira = list(.mor = "mor", fast = TRUE, center=FALSE),
  viper = list(verbose = FALSE,
               minsize = 0,
               .mor = "mor",
               .likelihood = "likelihood",
               pleiotropy = T,
               eset.filter = F),
  gsva = list(verbose = FALSE, method = "gsva"),
  ora = list(n_up=300, n_bottom=300, n_background=20000),
  fgsea = list(force_ties = T, options = list(nproc=4))
))

# Design
design <- tibble(
  set_name = 'dorothea', # name of the set resource
  bench_name = "dbd", # name of the benchmark data
  stats_list = stats_list,
  opts_list = opts_list,
  bexpr_loc = expr_fname, # benchmark data location
  bmeta_loc = meta_fname, # metadata location
  source_loc = netw_fname, # set source location
  source_col = "tf", # source name of the gene set source
  target_col = "target", # target name of the set source
  filter_col = "confidence", # column by which we wish to filter
  filter_crit = list(c('A','B','C')) # criteria by which we wish to filter
)

# Run benchmark
result <- run_benchmark(
  .design = design, # provide input tibble
  .minsize = 5, # filter gene sets with size < 10
  .form = TRUE, # format the benchmark results
  .perform = TRUE, # evaluate benchmarking performance
  .silent = FALSE, # silently run the pipeline
  .downsample_pr = TRUE, # downsample TNs for precision-recall curve
  .downsample_roc = TRUE, # downsample TNs for ROC
  .downsample_times = 100, # downsampling iterations (not used here)
  .url_bool = FALSE # whether to load from url
)

# Save result
saveRDS(result@bench_res, file.path(prc_path, 'rna_result.rds'))
