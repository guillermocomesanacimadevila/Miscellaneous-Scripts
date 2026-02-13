missingness_per_team <- function(number, path) {
  ID <- sprintf("DN%02d", number)
  file <- file.path(path, paste0(ID, ".csv"))
  df = read.csv(file)
  # missingness -> check for is.na OR NaN entries 
  na_count <- sum(is.na(df))
  return(na_count)
}

missingness_one_matrix <- function(in_dir, out_file) {
  files <- list.files(in_dir, full.names=TRUE, pattern="\\.csv$")
  feats <- c()
  for (f in files) {
    feats <- union(feats, names(read.csv(f, nrows=1)))
  }
  feats <- setdiff(feats, c("patient_id","team","X"))
  all <- data.frame()
  for (f in files) {
    d <- read.csv(f, stringsAsFactors=FALSE)
    X <- d[, feats, drop=FALSE]   
    miss <- is.na(X) | X == ""
    miss <- miss * 1
    out <- data.frame(patient_id=d$patient_id, team=d$team, miss, check.names=FALSE)
    all <- rbind(all, out)
  }
  write.csv(all, out_file, row.names=FALSE)
}
