#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(susieR) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript susie.R <pheno_id> <locus_sumstats> <L> <R_path> <n> <max_iter> [out_dir]")
}

pheno_id       <- args[1]
locus_sumstats <- args[2]
L              <- as.numeric(args[3])
R_path         <- args[4]
n              <- as.numeric(args[5])
max_iter       <- as.numeric(args[6])

locus_dir <- basename(dirname(locus_sumstats))
out_dir <- file.path("..", "outputs", "susie", pheno_id, locus_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
df <- read.table(locus_sumstats, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
R <- as.matrix(read.table(gzfile(R_path), header = FALSE))
storage.mode(R) <- "double"

run_susie <- function(df, R, L=1, n, estimate_residual_variance=FALSE, max_iter=5000) {
  z <- df$BETA / df$SE
  n_use <- median(df$N, na.rm = TRUE)
  
  # penalise R
  R <- (R + t(R))/2
  diag(R) <- 1
  lambda <- 1e-3
  R <- (1 - lambda) * R + lambda * diag(ncol(R))
  
  # cap L by effective rank of R
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  rankR <- sum(ev > 1e-4)
  L_use <- min(L, max(1, rankR - 1))
  message(sprintf("rankR=%d | requested L=%d | using L=%d", rankR, L, L_use))
  
  # susie
  fit <- susie_rss(
    z = z,
    R = R,
    n = n_use,
    L = L_use,
    estimate_residual_variance = estimate_residual_variance,
    max_iter = max_iter
  )
  
  cs_list <- fit$sets$cs
  cs_df <- data.frame()
  if (!is.null(cs_list) && length(cs_list) > 0) {
    for (i in seq_along(cs_list)) {
      nm <- names(cs_list)[i]
      if (is.null(nm) || nm == "") nm <- paste0("L", i)
      idx <- cs_list[[i]]
      df_cs <- data.frame(
        CS  = nm,
        SNP = df$SNP[idx],
        BP  = df$BP[idx],
        P   = df$P[idx],
        PIP = fit$pip[idx],
        stringsAsFactors = FALSE
      )
      df_cs <- df_cs[order(df_cs$PIP, decreasing=TRUE), ]
      cs_df <- rbind(cs_df, df_cs)
    }
  }
  
  list(fit = fit, cs = cs_df)
}

res <- run_susie(
  df, 
  R,
  L = L,
  n = n,
  max_iter = max_iter
)

pip_df <- data.frame(
  SNP = df$SNP,
  BP = df$BP,
  PIP = res$fit$pip,
  stringsAsFactors = FALSE
)

pip_df <- pip_df[order(pip_df$PIP, decreasing=TRUE), ]
saveRDS(res$fit, file.path(out_dir, paste0("res_", pheno_id, ".rds")))
write.table(res$cs, file.path(out_dir, paste0("cs_", pheno_id, ".tsv")), sep="\t", row.names=FALSE, quote=FALSE)
write.table(pip_df, file.path(out_dir, paste0("pip_", pheno_id, ".tsv")), sep="\t", row.names=FALSE, quote=FALSE)
