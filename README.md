# decoupleR_manuscript
Code to reproduce the results from decoupleR's manuscript

## Abstract
Computational methods allow the extraction of mechanistic signatures from omics 
data based on prior knowledge resources, reducing the dimensionality of the data
for increased statistical power and better interpretability. Here, we present 
decoupleR, a Bioconductor package containing different statistical methods to 
extract these signatures within a unified framework. decoupleR allows the user 
to flexibly test any method with any resource. It incorporates methods that take
into account the sign and weight of network interactions. Using decoupleR, we 
evaluated the performance of contemporary methods on transcriptomic and 
phospho-proteomic perturbation experiments.

## Benchmark pipeline
In this manuscript we have built a separate R package to benchmark decoupleR,
called decoupleRBench. You can install it by running:
```
devtools::install_github('saezlab/decoupleRBench')
```
You can check the source code in: https://github.com/saezlab/decoupleRBench

## Data
To retrieve all the necessary data, please run:
```
Rscript R/process/get_data.R
```

## Analysis
To reproduce all the analyses performed in the manuscript, please run:
```
# Run perturbation benchmark
Rscript R/process/run_rna_bench.R
Rscript R/process/run_php_bench.R
# Run noise analysis (This might take a while)
Rscript R/process/run_rna_noise.R #(~6 hours)
Rscript R/process/run_php_noise.R #(~2 hours)
```

## Figures
To generate all the figures shown in the manuscript, please run:
```
Rscript R/plot/supp_fig_1.R
Rscript R/plot/supp_fig_2.R
Rscript R/plot/supp_fig_3.R
Rscript R/plot/supp_fig_4.R
Rscript R/plot/supp_fig_5.R
Rscript R/plot/supp_fig_6.R
```

