get_Lambert <- function(){
  # Download Lambert list of TFs
  url <- 'https://ars.els-cdn.com/content/image/1-s2.0-S0092867418301065-mmc2.xlsx'
  fname <- file.path('data/', 'lambert.csv')
  download.file(url, fname)
  Sys.sleep(1)

  # Select the ones annotated as TF
  df <- readxl::read_xlsx(fname, sheet=2, skip=1)
  df <- dplyr::select(df, Name, `...4`)
  df <- dplyr::filter(df, `...4`=='Yes')
  df <- dplyr::select(df, Name)

  # Write
  write.csv2(df, file.path('data/', 'lambert.csv'), row.names=F)
}


filter_by_Lambert <- function(df, fname=file.path('data/', 'lambert.csv')){
  # Read TFs
  lambert <- read.csv2(fname)$Name
  dplyr::filter(df, tf %in% lambert)
}


get_regnetwork <- function(confs=c('High', 'Medium', 'Low')){
  # Download Regnetwork
  url <- 'http://www.regnetworkweb.org/export.jsp?format=csv&sql=SELECT+*+FROM+human+WHERE+%271%27+AND+%28confidence+%3D+%27{conf}%27%29+ORDER+BY+UPPER%28regulator_symbol%29+ASC'
  for (conf in confs) {
    get <- stringr::str_glue(url)
    fname <- file.path('data/', stringr::str_glue('regnetwork_{conf}.csv'))
    download.file(get, fname)
    Sys.sleep(1)
  }

  # Extract and merge
  df <- lapply(confs, function(conf){
    fname <- file.path('data/', stringr::str_glue('regnetwork_{conf}.csv'))
    df <- read.csv(fname)
    df
  })

  # Format
  df <- do.call(rbind, df)[,c('regulator_symbol','target_symbol','confidence')]
  colnames(df) <- c('tf', 'target', 'confidence')
  df['mor'] <- 1
  df['likelihood'] <- 1

  # Filter by Lambert
  df <- filter_by_Lambert(df)

  # Remove duplicates
  df <- df %>% distinct(tf, target, confidence, .keep_all = TRUE)

  # Write
  saveRDS(df, file.path('data/', 'regnetwork.rds'))

  # Remove tmp files
  for (conf in confs) {
    fname <- file.path('data/', stringr::str_glue('regnetwork_{conf}.csv'))
    file.remove(fname)
  }
}


unify_chea3_names <- function(tf_list) {
  # Insert '-' for NKX entries
  nkx_last_number <- tf_list %>% str_subset("NKX") %>% str_sub(-1)
  nkx_rest <- tf_list %>% str_subset("NKX") %>% str_sub(1,-2)
  nkx_corrected <- paste(nkx_rest, nkx_last_number, sep="-")
  # Correct NKX with '-' and mutate aliases
  tf_list %>%
    str_subset("NKX.*", negate = TRUE) %>%
    append(nkx_corrected) %>%
    str_replace("ZNF875", "HKR1") %>%
    str_replace("TBXT", "T") %>%
    str_replace("CBLL2","ZNF645") %>%
    str_replace("ZNF788P", "ZNF788") %>%
    str_replace("ZUP1", "ZUFSP")
}


get_chea3 <- function(){
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
    fname <- file.path('data/', stringr::str_glue('chea3_{name}.csv'))
    download.file(get, fname)
    Sys.sleep(1)
  }

  # Extract and merge
  df <- lapply(names, function(name){
    # Read file
    fname <- file.path('data/', stringr::str_glue('chea3_{name}.csv'))
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
  df <- filter_by_Lambert(df)

  # Remove duplicates
  df <- df %>% distinct(tf, target, confidence, .keep_all = TRUE)

  # Save
  saveRDS(df, file.path('data/', 'chea3.rds'))

  # Remove tmp files
  for (name in names) {
    fname <- file.path('data/', stringr::str_glue('chea3_{name}.csv'))
    file.remove(fname)
  }
}


get_dorothea <- function(){
  # Get dorothea
  data(dorothea_hs, package = "dorothea")

  # Filter by Lambert
  dorothea_hs <- filter_by_Lambert(dorothea_hs)

  # Remove duplicates
  dorothea_hs <- dorothea_hs %>% distinct(tf, target, confidence,
                                          .keep_all = TRUE)

  # Add likelihood
  dorothea_hs['likelihood'] <- 1

  # Save
  saveRDS(dorothea_hs, file.path('data/', 'dorothea.rds'))
}


get_data <- function(){
  dir.create('data', showWarnings = F)
  get_Lambert()
  get_dorothea()
  get_chea3()
  get_regnetwork()
}



