get_Lambert <- function(path){
  # Download Lambert list of TFs
  url <- 'https://ars.els-cdn.com/content/image/1-s2.0-S0092867418301065-mmc2.xlsx'
  fname <- file.path(path, 'lambert.csv')
  download.file(url, fname)
  Sys.sleep(1)

  # Select the ones annotated as TF
  df <- readxl::read_xlsx(fname, sheet=2, skip=1)
  df <- dplyr::select(df, Name, `...4`)
  df <- dplyr::filter(df, `...4`=='Yes')
  df <- dplyr::select(df, Name)

  # Write
  write.csv2(df, file.path(path, 'lambert.csv'), row.names=F)
}

filter_by_Lambert <- function(df, path){
  # Read TFs
  fname=file.path(path, 'lambert.csv')
  lambert <- read.csv2(fname)$Name
  dplyr::filter(df, tf %in% lambert)
}

get_dorothea <- function(path){
  # Get dorothea
  data(dorothea_hs, package = "dorothea")

  # Filter by Lambert
  dorothea_hs <- filter_by_Lambert(dorothea_hs, path)

  # Remove duplicates
  dorothea_hs <- dplyr::distinct(dorothea_hs, tf, target, confidence, .keep_all = TRUE)

  # Add likelihood
  dorothea_hs['likelihood'] <- 1

  # Filter by confidence
  dorothea_hs <- dplyr::filter(dorothea_hs, confidence %in% c("A","B","C"))

  # Save
  saveRDS(dorothea_hs, file.path(path, 'dorothea.rds'))
}


get_KSN <- function(path){
  # Download regulons
  url <- 'https://zenodo.org/record/5645208/files/kinase_regulons.txt?download=1'
  fname <- file.path(path, 'kinase_regulons.csv')
  download.file(url, fname)
  Sys.sleep(1)

  # Format
  KSN <- read.delim(fname)
  KSN <- tibble::tibble(source = KSN$Regulator,
                        likelihood = KSN$SS,
                        target = KSN$Target,
                        mor = 1,
                        confidence = "A")

  # Save
  saveRDS(KSN, file.path(path, 'KSN.rds'))

  # Remove tmp files
  file.remove(fname)
}

get_rna_data <- function(path){
  # Download
  expr_url <- 'https://zenodo.org/record/5645208/files/rna_expr.rds?download=1'
  meta_url <- 'https://zenodo.org/record/5645208/files/rna_meta.rds?download=1'
  download.file(expr_url, file.path(path, 'rna_expr.rds'))
  Sys.sleep(1)
  download.file(meta_url, file.path(path, 'rna_meta.rds'))
  Sys.sleep(1)

  # Filter by overlap of TFs
  rna_expr <- readRDS(file.path(path, 'rna_expr.rds'))
  rna_meta <- readRDS(file.path(path, 'rna_meta.rds'))
  network <- readRDS(file.path(path, 'dorothea.rds'))
  tfs <- network$tf
  rna_meta <- dplyr::filter(rna_meta, target %in% tfs)
  rna_expr <- rna_expr[, colnames(rna_expr) %in% rna_meta$id]

  # Save
  saveRDS(rna_expr, file.path(path, 'rna_expr.rds'))
  saveRDS(rna_meta, file.path(path, 'rna_meta.rds'))
}

get_php_data <- function(path){
  # Download
  expr_url <- 'https://zenodo.org/record/5645208/files/annotations.xlsx?download=1'
  meta_url <- 'https://zenodo.org/record/5645208/files/benchmark_data.xlsx?download=1'
  expr_fname <- file.path(path, 'php_expr.xlsx')
  meta_fname <- file.path(path, 'php_meta.xlsx')
  download.file(expr_url, expr_fname)
  Sys.sleep(1)
  download.file(meta_url, meta_fname)
  Sys.sleep(1)

  # Format expr
  phosphosites <- readxl::read_excel(expr_fname, sheet = "Phosphosites")
  logFC <- readxl::read_excel(expr_fname, sheet = "Foldchanges", na = "NA")
  benchmark_data <- tibble::tibble(
    cbind(
      ID = paste(phosphosites$ensp, phosphosites$positions, "P", sep = "_"), logFC)
    )
  benchmark_data <- tibble::column_to_rownames(benchmark_data, var = "ID")
  benchmark_data[is.na(benchmark_data)] <- 0

  # Filter out benchmamrk kinases not part of the KSN or that have target > 1
  metaData <- readxl::read_excel(meta_fname, sheet = "KinaseConditionPairs")
  KSN <- readRDS(file.path(path, 'KSN.rds'))
  tble <- table(metaData$Condition)
  ids <- names(tble)[tble == 1]
  benchmarkPKN <- metaData$Condition[(metaData$Kinase %in% KSN$source) & (metaData$Condition %in% ids)]
  benchmark_data <- benchmark_data[,colnames(benchmark_data) %in% benchmarkPKN]

  # Format meta
  metaData <- dplyr::rename(tibble::as_tibble(metaData),
                            id = Condition, target = Kinase, sign = Regulation)
  metaData <- dplyr::filter(metaData, id %in% colnames(benchmark_data))
  metaData <- dplyr::mutate(metaData, sign = replace(sign, sign == "up", 1))
  metaData <- dplyr::mutate(metaData, sign = replace(sign, sign == "down", -1))
  metaData$sign <- as.double(metaData$sign)

  # Save
  saveRDS(benchmark_data, file.path(path, 'php_expr.rds'))
  saveRDS(metaData, file.path(path, 'php_meta.rds'))

  # Remove tmp files
  file.remove(expr_fname)
  file.remove(meta_fname)
}

get_data <- function(path){
  # Create dir
  dir.create(path, showWarnings = F, recursive = T)

  # Get networks
  get_Lambert(path)
  get_dorothea(path)
  get_KSN(path)

  # Get benchmark data
  get_rna_data(path)
  get_php_data(path)
}

# Run
get_data(file.path('data','raw'))
