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


get_regnetwork <- function(path, confs=c('High', 'Medium', 'Low')){
  # Download Regnetwork
  url <- 'http://www.regnetworkweb.org/export.jsp?format=csv&sql=SELECT+*+FROM+human+WHERE+%271%27+AND+%28confidence+%3D+%27{conf}%27%29+ORDER+BY+UPPER%28regulator_symbol%29+ASC'
  for (conf in confs) {
    get <- stringr::str_glue(url)
    fname <- file.path(path, stringr::str_glue('regnetwork_{conf}.csv'))
    download.file(get, fname)
    Sys.sleep(1)
  }

  # Extract and merge
  df <- lapply(confs, function(conf){
    fname <- file.path(path, stringr::str_glue('regnetwork_{conf}.csv'))
    df <- read.csv(fname)
    df
  })

  # Format
  df <- do.call(rbind, df)[,c('regulator_symbol','target_symbol','confidence')]
  colnames(df) <- c('tf', 'target', 'confidence')
  df['mor'] <- 1
  df['likelihood'] <- 1

  # Filter by Lambert
  df <- filter_by_Lambert(df, path)

  # Remove miRNA targets
  df <- dplyr::filter(df, !grepl('hsa-', target))
  df <- dplyr::filter(df, !grepl('MIR', target))

  # Organize into layers:
  # High = High
  # Medium = High + Medium
  # Low = High + Medium + Low
  high_df <- dplyr::filter(df, confidence == 'High')

  medium_df <- dplyr::bind_rows(high_df, dplyr::filter(df, confidence == 'Medium'))
  medium_df$confidence <- 'Medium'

  low_df <- dplyr::bind_rows(medium_df, dplyr::filter(df, confidence == 'Low'))
  low_df$confidence <- 'Low'

  df <- dplyr::bind_rows(high_df, medium_df, low_df)
  df <- dplyr::distinct(df, tf, target, confidence, .keep_all = TRUE)

  # Write
  saveRDS(df, file.path(path, 'regnetwork.rds'))

  # Remove tmp files
  for (conf in confs) {
    fname <- file.path(path, stringr::str_glue('regnetwork_{conf}.csv'))
    file.remove(fname)
  }
}


unify_chea3_names <- function(tf_list) {
  # Insert '-' for NKX entries
  nkx_last_number <- stringr::str_sub(stringr::str_subset(tf_list, "NKX"), -1)
  nkx_rest <- stringr::str_sub(stringr::str_subset(tf_list, "NKX"), 1,-2)
  nkx_corrected <- paste(nkx_rest, nkx_last_number, sep="-")
  # Correct NKX with '-' and mutate aliases
  tf_list <- stringr::str_subset(tf_list, "NKX.*", negate = TRUE)
  tf_list <- append(tf_list, nkx_corrected)
  tf_list <- stringr::str_replace(tf_list, "ZNF875", "HKR1")
  tf_list <- stringr::str_replace(tf_list, "TBXT", "T")
  tf_list <- stringr::str_replace(tf_list, "CBLL2","ZNF645")
  tf_list <- stringr::str_replace(tf_list, "ZNF788P", "ZNF788")
  tf_list <- stringr::str_replace(tf_list, "ZUP1", "ZUFSP")
  tf_list
}


get_chea3 <- function(path){
  # Download chea3
  url <- 'https://maayanlab.cloud/chea3/assets/tflibs/{name}.gmt'
  names <- c(
    'ARCHS4_Coexpression',
    'ENCODE_ChIP-seq',
    'Enrichr_Queries',
    'GTEx_Coexpression',
    'Literature_ChIP-seq',
    'ReMap_ChIP-seq'
  )
  for (name in names) {
    get <- stringr::str_glue(url)
    fname <- file.path(path, stringr::str_glue('chea3_{name}.csv'))
    download.file(get, fname)
    Sys.sleep(1)
  }

  # Extract and merge
  df <- lapply(names, function(name){
    # Read file
    fname <- file.path(path, stringr::str_glue('chea3_{name}.csv'))
    f <- file(fname,open="r")
    lines <-readLines(f)
    close(f)

    dfs <- lapply(lines, function(line){
      # Remove tabs
      line <- stringr::str_split(line, '\t')[[1]]

      # Select first element as TF
      tf <- stringr::str_split(line[1], '_')[[1]][1]

      # The rest are targets
      target <- line[2:length(line)]

      # Repeat elements
      tf <- rep(tf, length(target))
      confidence <- rep(name, length(target))

      # Save as df
      data.frame('tf'=tf, 'confidence'=confidence, 'target'=target)
    })
    dfs <- do.call(rbind, dfs)
  })

  # Format
  df <- do.call(rbind, df)
  df['mor'] <- 1
  df['likelihood'] <- 1

  # Change mislabeled TFs
  df[['tf']] <- unify_chea3_names(df[['tf']])

  # Filter by Lambert
  df <- filter_by_Lambert(df, path)

  # Remove duplicates
  df <- dplyr::distinct(df, tf, target, confidence, .keep_all = TRUE)

  # Save
  saveRDS(df, file.path(path, 'chea3.rds'))

  # Remove tmp files
  for (name in names) {
    fname <- file.path(path, stringr::str_glue('chea3_{name}.csv'))
    file.remove(fname)
  }
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

  # Save
  saveRDS(dorothea_hs, file.path(path, 'dorothea.rds'))
}


get_KSN <- function(path){
  # Download regulons
  url <- 'https://docs.google.com/uc?export=download&id=1vkSIuIneXb30VrMuFzGd2yWofQl962cT'
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
  expr_url <- 'https://zenodo.org/record/4322914/files/dorothea_bench_expr.rds?download=1'
  meta_url <- 'https://zenodo.org/record/4322914/files/dorothea_bench_meta.rds?download=1'
  download.file(expr_url, file.path(path, 'rna_expr.rds'))
  Sys.sleep(1)
  download.file(meta_url, file.path(path, 'rna_meta.rds'))
  Sys.sleep(1)
}

get_php_data <- function(path){
  # Download
  expr_url <- 'https://docs.google.com/uc?export=download&id=1XtnPk5VhrJ5_Vn9J0Roo5rbkTNYLsw4G'
  meta_url <- 'https://docs.google.com/uc?export=download&id=1N8nsCUuxRZEUC6sIssyyROBk7FwAZfc_'
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

  # Filter out benchmamrk kinases not part of the KSN
  metaData <- readxl::read_excel(meta_fname, sheet = "KinaseConditionPairs")
  KSN <- readRDS(file.path(path, 'KSN.rds'))
  benchmarkPKN <- unique(metaData$Condition[metaData$Kinase %in% KSN$source])
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
  get_chea3(path)
  get_regnetwork(path)
  get_KSN(path)

  # Get benchmark data
  get_rna_data(path)
  get_php_data(path)
}


# Run
get_data('data/raw')

