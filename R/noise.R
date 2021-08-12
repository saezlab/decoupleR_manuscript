library(decoupleRBench)
library(dplyr)
library(purrr)

bexpr_loc = file.path('data', "dorothea_bench_expr.rds")
bmeta_loc <- file.path('data', "dorothea_bench_meta.rds")
# List of the methods to call
stats_list = list(c("mean", "pscira","viper"))
# List of options for each method
opts_list <- list(list(
  #ora = list(n_up=300, n_bottom=300, n_background=20000),
  mean = list(.mor = "mor"),
  pscira = list(.mor = "mor"),
  viper = list(verbose = FALSE,
               minsize = 0,
               .mor = "mor",
               .likelihood = "likelihood"
               )
))


source_loc <- file.path('data','dorothea.rds')
design <- tibble(
  set_name = c('dorothea','dorothea'),#, 'dorothea_add50%', 'dorothea_del50'),
  bench_name = c("dbd",'dbd'),#,"dbd","dbd"),
  stats_list = c(stats_list,stats_list),#,stats_list,stats_list),
  opts_list = c(opts_list,opts_list),#,opts_list,opts_list),
  bexpr_loc = c(bexpr_loc,bexpr_loc),#,bexpr_loc,bexpr_loc),
  bmeta_loc = c(bmeta_loc,bmeta_loc),#,bmeta_loc,bmeta_loc),
  source_loc = c(source_loc,source_loc),#,source_loc,source_loc),
  source_col = c("tf",'tf'),#,"tf","tf"),
  target_col = c("target",'target'),#,"target","target"),
  filter_col = c("confidence",'confidence'),#,"confidence","confidence"),
  filter_crit = list(c("A"),c("A", 'B', 'C')),#,c("A", 'B', 'C')),
  noise_crit = list(NA,
                    NA
                    #list(mode='add', perc=0.5, seed=1996),
                    #list(mode='del', perc=0.5, seed=1996)
                    )
)

estimate <- run_benchmark(
  .design = design, # provide input tibble
  .minsize = 5, # filter gene sets with size < 10
  .form = TRUE, # format the benchmark results
  .perform = TRUE, # evaluate benchmarking performance
  .silent = FALSE, # silently run the pipeline
  .downsample_pr = FALSE, # downsample TNs for precision-recall curve
  .downsample_roc = TRUE, # downsample TNs for ROC
  .downsample_times = 100, # downsampling iterations (not used here)
  .url_bool = FALSE # whether to load from url
)


estimate@bench_res$activity %>%
  map(function(df){

  })

# Get percentage of correct directions
estimate@bench_res$activity[[3]] %>%
  filter(tf == target) %>%
  mutate(same_sign = sign(score) == sign(sign)) %>%
  summarise(perc = sum(same_sign) / n())

result <- estimate@summary$summary_table %>%
  mutate('type' = 'estimate') %>%
  select(type, statistic, auc)






