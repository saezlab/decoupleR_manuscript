library(decoupleR)
suppressPackageStartupMessages(library(tidyverse))
source(file.path('R', 'process', 'methods_params.R'))

options( warn = -1 )

# Parse args
args <- commandArgs(trailingOnly = TRUE)
mat_path <- args[[1]]
net_path <- args[[2]]
method <- args[[3]]

if (method == 'gsea'){
  method <- 'fgsea'
}

# Methods dictionary
methods <- list(
  aucell = run_aucell,
  udt = run_udt,
  mdt = run_mdt,
  wmean = run_wmean,
  ulm = run_ulm,
  mlm = run_mlm,
  wsum = run_wsum,
  viper = run_viper,
  gsva = run_gsva,
  ora = run_ora,
  fgsea = run_fgsea
)

# Open
mat <- as.matrix(read.csv(mat_path))
colnames(mat) <- as.character(1:ncol(mat))
rownames(mat) <- as.character(1:nrow(mat))
net <- as_tibble(read.csv(net_path)) %>%
  mutate(source=as.character(source), target=as.character(target))

# Run
f <- methods[[method]]
params <- opts_list[[1]][[method]]
params[['mat']] <- mat
params[['net']] <- net
res <- do.call(f, params)
