library(decoupleRBench)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)


bexpr_loc = file.path('data', "dorothea_bench_expr.rds")
# benchmark metadata
bmeta_loc <- file.path('data', "dorothea_bench_meta.rds")

# List of the methods to call
stats_list = list(c('mean','pscira','scira','viper','gsva','ora','fgsea'))
stats_list = list(c('mean', 'viper'))

# List of options for each method
opts_list <- list(list(
  mean = list(times=100, .mor = "mor"),
  pscira = list(times=100, .mor = "mor"),
  scira = list(.mor = "mor", fast = FALSE, center=FALSE),
  viper = list(verbose = FALSE,
               minsize = 0,
               .mor = "mor",
               .likelihood = "likelihood"),
  gsva = list(verbose = FALSE, method = "gsva"),
  ora = list(n_up=300, n_bottom=300, n_background=20000),
  fgsea = list(force_ties = T, options = list(nproc=4))
))
opts_list <- list(list(
  mean = list(times=10, .mor = "mor"),
  viper = list(verbose = FALSE,
               minsize = 0,
               .mor = "mor",
               .likelihood = "likelihood")
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
input_tibble <- input_tibble %>%
  filter(set_name == 'dorothea')

input_tibble <- input_tibble['dorothea.3',]

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
saveRDS(dor_run, 'out.rds')

df <- dor_run@bench_res %>% dplyr::select(set_name, filter_crit, activity) %>% tidyr::unnest(activity)
df <- df %>% pivot_wider(names_from=statistic, values_from=score)

df <- dor_run@bench_res %>%
  unnest(roc) %>%
  select(.data$set_name, .data$bench_name, .data$filter_crit,
         .data$statistic, .data$auc, .data$coverage) %>%
  distinct() %>%
  unite("name_lvl", .data$set_name, .data$bench_name, .data$filter_crit)


ggplot(df, aes(x=name_lvl , y=statistic, color=auc, size=coverage)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_colour_gradient2(low = 'blue', high='red', mid='gray', midpoint = 0.5)


# Correlation
df <- out@bench_res %>%
  dplyr::select(set_name, filter_crit, activity) %>%
  tidyr::unnest(activity) %>%
  select(set_name, filter_crit, id, tf, statistic, score) %>%
  pivot_wider(names_from=statistic, values_from=score) %>%
  select(-set_name, -filter_crit, -id, -tf)

corr_matrix <- matrix(0, ncol(df), ncol(df))
colnames(corr_matrix) <- colnames(df)
rownames(corr_matrix) <- colnames(df)

for (name_a in colnames(df)) {
  for (name_b in colnames(df)) {
    corr_matrix[name_a,name_b] <- cor(abs(df[[name_a]]), abs(df[[name_b]]))
  }
}
pheatmap(corr_matrix, display_numbers=T)

# Jaccard indexes
df <- out@bench_res %>%
  mutate(set_name = paste0(set_name, '_', filter_crit)) %>%
  select(set_name, activity) %>%
  tidyr::unnest(activity) %>%
  select(set_name, id, tf, statistic, score) %>%
  mutate(score = abs(score)) %>%
  arrange(set_name, id, statistic, desc(score)) %>%
  group_by(set_name, id, statistic) %>%
  slice_head(n=1000)


x <- df %>%
  group_by(set_name, id, statistic) %>%
  summarize(tf = list(tf), .groups='drop')

statistics <- unique(df$statistic)

vec <- c()
for (i in 5) {
  for (j in 5) {
    vec[i,j] <- 1
  }
}


x %>%
  filter(set_name == 'dorothea_A') %>%
  group_by(set_name, id) %>%
  group_split() %>%
  map(function(df){
    n_stats <- nrow(df)
    for (i in n_stats) {
      stat_a <- df[i,'statistic']
      for (j in n_stats) {
        if (i < j) {
          stat_b <- df[j,'statistic']
          df %>%
            summarize(Jaccard = jaccard_index())
        }
      }
    }
  })




statistics <- unique(df$statistic)
networks <- unique(df$set_name)
samples <- unique(df$id)

pairs <- c()
l <- 1
for (i in 1:length(statistics)) {
  for (j in 1:length(statistics)) {
    if (j < i) {
      pairs[l] <- paste0(statistics[i],'_',statistics[j])
      l = l + 1
    }
  }
}

n_tops <- 1:10 * 10

jaccard_index <- function(a,b){
  return(length(intersect(a, b)) / length(union(a, b)))
}


tmp <- purrr::map(networks, function(net){
  purrr::map(samples, function(smp){
    purrr::map(pairs, function(pair){
      pair <- stringr::str_split(pair, '_')
      stat_a <- pair[1]
      stat_b <- pair[2]
      a <- df %>%
        filter(set_name == net, id == smp, statistic == stat_a) %>%
        slice_head(n=10)
      b <- df %>%
        filter(set_name == net, id == smp, statistic == stat_b) %>%
        slice_head(n=10)
      jacc <- jaccard_index(a$tf, b$tf)
      tibble(
        n_top = 10,
        network = net,
        sample = smp,
        stat_a = stat_a,
        stat_b = stat_b,
        jacc = jacc
      )
    })
  })
})







df %>%
  group_by(set_name) %>%
  group_split() %>%
  purrr:::map(function(net_df){
    net_df %>%
    group_by(id) %>%
    group_split() %>%
    purrr::map(function(smp_df){
      smp_df %>%
        group_by(statistic) %>%
        slice_head(n=10) %>%
        purrr::map(function(){

      })
    })
  })
  slice_head(n=10) %>%
  summarise(Union=length(union(tf)))




