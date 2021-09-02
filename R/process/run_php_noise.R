library(decoupleRBench)
library(dplyr)
library(purrr)
library(tidyverse)

# Paths to benchmark data, benchmark metadata and kinase substrate network
# Benchmark data contains 82 perturbation experiments covering 27 unique kinases
# Network contains 92 kinases with regulon size of 57 ± 86
raw_path <- file.path('data', 'raw')
prc_path <- file.path('data', 'prc')
dir.create(prc_path, showWarnings = F, recursive = T)
expr_fname <- file.path(raw_path, "php_expr.rds")
meta_fname <- file.path(raw_path, "php_meta.rds")
netw_fname <- file.path(raw_path, 'KSN.rds')

# List of the methods to call
stats_list = list(c('mean','pscira','scira','viper','gsva','ora','fgsea'))

# List of options for each method
opts_list <- list(list(
  mean = list(times=100, .mor = "mor"),
  pscira = list(times=100, .mor = "mor"),
  scira = list(.mor = "mor", fast = FALSE, center=FALSE),
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
  set_name = 'normal', # name of the set resource
  bench_name = '0', # name of the benchmark data
  stats_list = stats_list,
  opts_list = opts_list,
  bexpr_loc = expr_fname, # benchmark data location
  bmeta_loc = meta_fname, # metadata location
  source_loc = netw_fname, # set source location
  source_col = "source", # source name of the gene set source
  target_col = "target", # target name of the set source
  filter_col = "confidence", # column by which we wish to filter
  filter_crit = list(c('A')), # criteria by which we wish to filter,
  noise_crit = list(NA)
)

design <- bind_rows(
  design,
  map_df(c('add', 'del'), function(mode){
    map_df(c(0.25, 0.5, 0.75), function(perc){
      map_df(1:5, function(i){
        design %>%
          mutate(set_name=paste0(mode,as.character(i)), bench_name=as.character(perc),
                 noise_crit = list(list(mode=mode, perc=perc, seed=i)))
      })
    })
  })
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
  .downsample_times = 5, # downsampling iterations
  .url_bool = FALSE # whether to load from url
)

# Save result
saveRDS(result@bench_res, file.path(prc_path, 'php_noise.rds'))
