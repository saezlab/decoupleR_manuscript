library(decoupleRBench)
library(dplyr)
library(purrr)
library(tidyverse)

# Paths to benchmark data, benchmark metadata and kinase substrate network
# Benchmark data contains 82 perturbation experiments covering 27 unique kinases
# Network contains 92 kinases with regulon size of 57 ± 86
bexample_url <- file.path('data', "benchmark_data.rds")
bmeta_url <- file.path('data', "benchmark_metaData.rds")
source_url <- file.path('data', "kinase_substrate_network.rds")

# Design contains statistical methods that take weights into account
design_row <-
  tibble(
    set_name = "KSN_weighted", # name of the set resource
    bench_name = "beltrao", # name of the benchmark data
    stats_list = list( # a list of the stats to call
      c(
        "scira",
        "pscira",
        "mean",
        "viper"
      )
    ),
    opts_list = list(list( # list of options for each stat method
      scira = list(.mor = "mor",
                   .likelihood = "likelihood"),
      pscira = list(times = 100,
                    .mor = "mor",
                    .likelihood = "likelihood"),
      mean = list(times = 100,
                  .mor = "mor",
                  .likelihood = "likelihood"),
      viper = list(verbose = FALSE,
                   minsize = 0,
                   .mor = "mor",
                   .likelihood = "likelihood")
    )),
    bexpr_loc = bexample_url, # benchmark data location
    bmeta_loc = bmeta_url, # metadata location
    source_loc = source_url, # set source location
    source_col = "source", # source name of the gene set source
    target_col = "target", # target name of the set source'
    filter_col = "confidence", # column by which we wish to filter
    filter_crit = list(c("A")) # criteria by which we wish to filter
  )

# input tibble for one run with the weighted network
# and one where the weight is removed
input_tibble <- bind_rows(
  design_row,
  design_row %>%
    mutate(set_name ="KSN_unweighted") %>%
    mutate(weight_crit = list(list(.likelihood = "likelihood")))
)


# run decoupleRBenchmark
estimate <- run_benchmark(
  .design = input_tibble, # provide input tibble
  .minsize = 5, # filter gene sets with size < 10
  .form = TRUE, # format the benchmark results
  .perform = TRUE, # evaluate benchmarking performance
  .silent = TRUE, # silently run the pipeline
  .downsample_pr = FALSE, # downsample TNs for precision-recall curve
  .downsample_roc = TRUE, # downsample TNs for ROC
  .downsample_times = 100, # downsampling iterations
  .url_bool = FALSE # whether to load from url
)


# extract each auc per permutation run
auc_downsampling <- lapply(estimate@bench_res$roc,
                           function(x) x %>% group_by(run) %>% summarize(raw_auc = unique(raw_auc)) %>% pull(raw_auc))

estimate@bench_res <- add_column(estimate@bench_res, auc_downsampling)

# get indices of comparison groups (weighted and unweighted network)
comp_groups <- estimate@bench_res %>% group_by(statistic) %>% group_rows()
names(comp_groups) <- estimate@bench_res %>% pull(statistic) %>% unique()

# perform one tailed t-test between comparison groups
t.test_res <- map(comp_groups, function(group_idx) {
  t.test(estimate@bench_res$auc_downsampling[group_idx][[1]],
         estimate@bench_res$auc_downsampling[group_idx][[2]],
         alternative="greater")
})

# result table containing t- and p-values for each method comparing weighted and unweighted networks
results <- tibble(statistic = names(t.test_res),
                  t.value = unlist(lapply(t.test_res, `[`, "statistic")),
                  p.value = unlist(lapply(t.test_res, `[`, "p.value"))) %>% arrange(p.value)
show(results)

# boxplot sorted by methods with lowest p-value (t-test)
boxplot_tibble <- bind_cols(auc = unlist(auc_downsampling),
                            statistic = rep(estimate@bench_res$statistic, each = 100),
                            network = rep(estimate@bench_res$set_name, each = 100))
boxplot_tibble$statistic <- factor(boxplot_tibble$statistic, levels = results$statistic)
boxplot_tibble$network <- factor(boxplot_tibble$network, levels = c("KSN_weighted", "KSN_unweighted"))

ggplot(boxplot_tibble,aes(fill = network, x = statistic, y = auc)) +
    geom_boxplot()

