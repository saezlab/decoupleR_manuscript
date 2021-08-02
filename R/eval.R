library(decoupleRBench)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)


bexpr_loc = file.path('data', "dorothea_bench_expr.rds")
# benchmark metadata
bmeta_loc <- file.path('data', "dorothea_bench_meta.rds")

# List of the methods to call
stats_list = list(c("mean","pscira","scira","viper","gsva","ora"))
stats_list = list(c("pscira","ora"))

# List of options for each method
opts_list <- list(list(
  scira = list(.mor = "mor"),
  pscira = list(.mor = "mor"),
  mean = list(.mor = "mor"),
  viper = list(verbose = FALSE,
               minsize = 0,
               .mor = "mor",
               .likelihood = "likelihood"),
  gsva = list(verbose = FALSE, method = "gsva"),
  ora = list()
))

opts_list <- list(list(
  pscira = list(.mor = "mor"),
  ora = list()
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
      filter_crit = filter_crit # criteria by which we wish to filter
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
  .downsample_pr = FALSE, # downsample TNs for precision-recall curve
  .downsample_roc = FALSE, # downsample TNs for ROC
  .downsample_times = 100, # downsampling iterations (not used here)
  .url_bool = FALSE # whether to load from url
)

design_row <- tibble(
    set_name = "dorothea", # name of the set resource
    bench_name = "dbd", # name of the benchmark data
    stats_list = list( # a list of the stats to call
      c(
        "pscira",
        "ora"
      )
    ),
    opts_list = list(list( # list of options for each stat method
      pscira = list(.mor = "mor"),
      ora = list()
    )),
    bexpr_loc = bexpr_loc, # benchmark data location
    bmeta_loc = bmeta_loc, # metadata location
    source_loc = file.path('data','dorothea.rds'), # set source location
    source_col = "tf", # source name of the gene set source
    target_col = "target", # target name of the set source
    filter_col = "confidence", # column by which we wish to filter
    filter_crit = list(c("A")) # criteria by which we wish to filter
  )

dor_run <- run_benchmark(
  .design = design_row, # provide input tibble
  .minsize = 10, # filter gene sets with size < 10
  .form = TRUE, # format the benchmark results
  .perform = TRUE, # evaluate benchmarking performance
  .silent = FALSE, # silently run the pipeline
  .downsample_pr = FALSE, # downsample TNs for precision-recall curve
  .downsample_roc = FALSE, # downsample TNs for ROC
  .downsample_times = 100, # downsampling iterations (not used here)
  .url_bool = FALSE # whether to load from url
)








library(decoupleR)
library(dplyr)
library(tibble)

inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
thr <- 0.01
n_background <- NULL
with_ties <- TRUE

regulons <- decoupleR::convert_to_ora(network, tf, target)

ora_check_ns <- function(threshold, n_background, network, mat) {
  if (is.null(n_background)) {
    n_background <- network %>%
      pull(.data$target) %>%
      unique() %>%
      union(rownames(mat)) %>%
      length()
  } else if (n_background < 0) {
    abort("`n` must be a non-missing positive number.")
  }

  n_up <- ceiling(threshold * length(mat))

  c(n_up, n_background)
}
ns <- ora_check_ns(thr, n_background, network, mat)
n_up <- ns[1]
n_background <- ns[2]
ora_slice_targets <- function(mat, n_up, with_ties) {
  mat %>%
    as_tibble(rownames = "target") %>%
    tidyr::pivot_longer(
      cols = -.data$target,
      names_to = "condition",
      values_to = "value"
    ) %>%
    group_by(.data$condition) %>%
    slice_max(., abs(.data$value), n = n_up, with_ties = with_ties) %>%
    summarise(
      targets = purrr::set_names(list(.data$target), .data$condition[1]),
      .groups = "drop"
    ) %>%
    pull(.data$targets)
}
targets <- ora_slice_targets(mat, n_up, with_ties)
expand.grid(tf = names(regulons), condition = names(targets)) %>%
  rowwise(.data$tf, .data$condition)


run_hyper_test <- function(top_genes, target_genes, n_background){
  GP <- length(target_genes)
  intersection <- intersect(top_genes, target_genes)
  GL <- length(intersection)
  p_value <- phyper(q=GL - 1, m=GP, n=n_background-GP, k=length(top_genes),
         lower.tail = FALSE, log.p = FALSE)
  p_value
}

run_hyper_test(targets$GSM2753335, regulons$FOXO4, n_background)
