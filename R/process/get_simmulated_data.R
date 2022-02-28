library(tidyverse)

# Parameters for simulated data
n_features <- 20000
n_targets <- 40
vals <- list('1'=10, '2'=50, '3'=250, '4'=1250, '5'=2500)
out_path <- file.path('data','raw')

withr::with_seed(42, {
  for (i in 1:length(vals)) {
    i <- as.character(i)

    # Create matrix
    n_samples <- vals[[i]]
    mat <- matrix(stats::rnorm(n_features*n_samples), nrow=n_features)
    mat[abs(mat) < 1] <- 0.0

    # Save
    write.csv(mat, file.path(out_path, paste0(i,'_mat.csv')), row.names = F)
  }

  # Create net
  i <- '3'
  net <- tibble(source=as.character(1:vals[[i]])) %>%
    mutate(target=map(source, function(i){
      as.character(sample.int(n_features, size = n_targets, replace = F))
    })) %>%
    unnest(target)
  net[['mor']] <- rnorm(nrow(net))

  # Save
  write.csv(net, file.path(out_path, paste0('net.csv')), row.names = F)
})



