library(decoupleRBench)
library(dplyr)
library(purrr)
library(tidyverse)
source(file.path('R', 'process', 'methods_params.R'))

# Paths to benchmark data, benchmark metadata and network
raw_path <- file.path('data', 'raw')
prc_path <- file.path('data', 'prc')
dir.create(prc_path, showWarnings = F, recursive = T)
expr_fname <- file.path(raw_path, "rna_expr.rds")
meta_fname <- file.path(raw_path, "rna_meta.rds")
netw_fname <- file.path(raw_path, 'dorothea.rds')

# Design
design <- tibble(
  set_name = 'normal', # name of the set resource
  bench_name = '0', # name of the benchmark data
  stats_list = stats_list,
  opts_list = opts_list,
  bexpr_loc = expr_fname, # benchmark data location
  bmeta_loc = meta_fname, # metadata location
  source_loc = netw_fname, # set source location
  source_col = "tf", # source name of the gene set source
  target_col = "target", # target name of the set source
  filter_col = "confidence", # column by which we wish to filter
  filter_crit = list(c('A','B','C')), # criteria by which we wish to filter,
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
  .downsample_times = 20, # downsampling iterations
  .url_bool = FALSE # whether to load from url
)

# Save result
result <- result@bench_res %>%
  mutate(
    activity = map(activity, function(df){
      df %>%
        select(run_id, statistic, source, id, score, statistic_time)
    }),
    roc = map(roc, function(df){
      df %>%
        distinct(run, .keep_all=T) %>%
        select(raw_auc, run, auc)
    }),
    prc = map(prc, function(df){
      df %>%
        distinct(run, .keep_all=T) %>%
        select(raw_auc, run, auc)
    })
  )

saveRDS(result, file.path(prc_path, 'rna_noise.rds'))
