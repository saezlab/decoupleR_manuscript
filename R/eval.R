library(decoupleRBench)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)


bexpr_loc = file.path('data', "dorothea_bench_expr.rds")
# benchmark metadata
bmeta_loc <- file.path('data', "dorothea_bench_meta.rds")

# List of the methods to call
stats_list = list(c("mean","pscira","scira","viper","gsva","ora","fgsea"))

# List of options for each method
opts_list <- list(list(
  mean = list(times=10, .mor = "mor"),
  pscira = list(times=10, .mor = "mor"),
  scira = list(.mor = "mor"),
  viper = list(verbose = FALSE,
               minsize = 0,
               .mor = "mor",
               .likelihood = "likelihood"),
  gsva = list(verbose = FALSE, method = "gsva"),
  ora = list(n_up=300, n_bottom=300, n_background=20000),
  fgsea = list(force_ties = T)
))

# List of path
net_paths <- list(
  'dorothea' = file.path('data', 'dorothea.rds'),
  'chea3'= file.path('data', 'chea3.rds'),
  'regnetwork'= file.path('data', 'regnetwork.rds')
  )

# Filtering criteria
net_confs <- list(
  'dorothea' = list(
    c('A'),
    c('A','B','C'),
    c('A','B','C','D','E')
    ),
  'chea3' = list(
    c('ARCHS4_Coexpression'),
    c('ENCODE_ChIP-seq'),
    c('Enrichr_Queries'),
    c('GTEx_Coexpression'),
    c('Literature_ChIP-seq'),
    c('ReMap_ChIP-seq')
  ),
  'regnetwork' = list(
    c('High'),
    c('Medium'),
    c('Low')
  )
)

# Generate input tibble
input_tibble <- lapply(names(net_confs), function(network){
  tmp <- lapply(net_confs[network], function(filter_crit){
    design_row <- tibble(
      set_name = network, # name of the set resource
      bench_name = "dbd", # name of the benchmark data
      stats_list = stats_list,
      opts_list = opts_list,
      bexpr_loc = bexpr_loc, # benchmark data location
      bmeta_loc = bmeta_loc, # metadata location
      source_loc = net_paths[[network]], # set source location
      source_col = "tf", # source name of the gene set source
      target_col = "target", # target name of the set source
      filter_col = "confidence", # column by which we wish to filter
      filter_crit = filter_crit, # criteria by which we wish to filter
      noise_crit = NA
    )
  })
  tmp <- do.call(rbind, tmp)
})
input_tibble <- do.call(rbind, input_tibble)

# Run benchmark
dor_run <- run_benchmark(
  .design = input_tibble, # provide input tibble
  .minsize = 10, # filter gene sets with size < 10
  .form = TRUE, # format the benchmark results
  .perform = TRUE, # evaluate benchmarking performance
  .silent = FALSE, # silently run the pipeline
  .downsample_pr = TRUE, # downsample TNs for precision-recall curve
  .downsample_roc = TRUE, # downsample TNs for ROC
  .downsample_times = 100, # downsampling iterations (not used here)
  .url_bool = FALSE # whether to load from url
)
