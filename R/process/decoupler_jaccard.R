library(tidyverse)

out <- readRDS("out.rds")
decoupler_out <- out@bench_res %>%
  mutate(set_name = paste0(set_name, '_', filter_crit)) %>%
  select(set_name, activity) %>%
  tidyr::unnest(activity) %>%
  select(set_name, id, tf, statistic, score) %>%
  mutate(score = abs(score)) %>%
  arrange(set_name, id, statistic, desc(score))

# Jaccard ix function -------------------------

get_pairwise_jaccard <- function(dat) {

  dat %>%
    mutate(j_ix_res = map(statistic, function(stat) {

      init_df <- dat %>%
        filter(statistic == stat)

      query <- dat %>%
        filter(statistic != stat)

      init_TF <- init_df[["tf"]][[1]]

      query <- query %>%
        mutate(j_ix = map_dbl(tf, function(q_tf) {

          length(intersect(q_tf, init_TF)) / length(union(q_tf, init_TF))

        })) %>%
        dplyr::select(statistic, j_ix)

    })) %>%
    dplyr::select(statistic, j_ix_res)

}


# Network size ------------------------------
network_size_meta <- decoupler_out %>%
  group_by(set_name) %>%
  summarize(net_size = length(unique(tf)))

# Get top n per percent
seq_perc <- c(0.05, 0.1, 0.15,0.2, 0.25, 0.3)
seq_perc <- set_names(seq_perc)

topn_list <- map(seq_perc, function(perc) {
  network_size_meta %>%
    mutate(top_n = (net_size * perc) %>% ceiling())
})


top_tf_list <- map(topn_list, function(top_df) {

  decoupler_out %>%
    group_by(set_name, id, statistic) %>%
    nest() %>%
    left_join(top_df, by = "set_name")

})

# This is your top n TF dataframe
library(furrr)
plan(multisession, workers = 8)

top_tfs <- furrr::future_map(top_tf_list, function(top_df){

  top_df %>%
    mutate(tf = map2(data, top_n, function(dat, n) {

      dat %>%
        dplyr::slice(1:n) %>%
        pull(tf) %>%
        as.character()

    })) %>%
    dplyr::select(set_name, id, statistic, tf)

})

# Jaccard ix
plan(multisession, workers = 8)
jaccard_list <- furrr::future_map(top_tfs, function(decoupler_out) {
  decoupler_out %>%
    group_by(set_name, id) %>%
    nest() %>%
    mutate(exp_net_jix = map(data, get_pairwise_jaccard)) %>%
    dplyr::select(-data) %>%
    unnest()
})

# FINAL PRODUCT ---------------

df <- jaccard_list %>% enframe() %>% unnest() %>% unnest()

df %>%
  mutate(name = as.double(name)) %>%
  filter(statistic == 'viper' | statistic == 'normalized_mean') %>%
  ggplot(aes(x=name, y=j_ix, color=set_name)) + geom_point() +geom_smooth(method='lm', aes(fill=set_name))

df %>%
  mutate(name = as.double(name)) %>%
  filter(statistic == 'scira' | statistic == 'mean') %>%
  ggplot(aes(x=name, y=j_ix, color=set_name)) + geom_point() +geom_smooth(method='lm', aes(fill=set_name))

df %>%
  mutate(name = as.double(name)) %>%
  filter(statistic == 'scira' | statistic == 'mean')

df %>%
  mutate(name = as.double(name)) %>%
  filter(statistic == 'mean' & statistic1 == 'pscira') %>%
  ggplot(aes(x=name, y=j_ix, color=set_name)) + geom_smooth(method='auto')

df %>%
  mutate(name = as.double(name)) %>%
  filter(statistic == 'viper' & statistic1 == 'normalized_mean') %>%
  ggplot(aes(x=name, y=j_ix, color=set_name)) + geom_smooth(method='auto')

df %>%
  mutate(name = as.double(name)) %>%
  ggplot(aes(x=name, y=j_ix, color=set_name)) + geom_smooth(method='glm') + theme_classic() +ylim(0,1)


df %>%
  mutate(name = as.double(name)) %>%
  ggplot(., aes(x = name, y = j_ix)) +
  geom_smooth(method="auto", se=TRUE) +
  labs(y = "Jaccard index", x = "% of top TFs")

