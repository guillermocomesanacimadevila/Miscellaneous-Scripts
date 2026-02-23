#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(susieR)
  library(Matrix)
})

# args 
# pheno_id
# R_matrix
# locus file 
# L param
# n_iter
# out_dir

run_susie <- function(df, R, L=1, n=NULL, estimate_residual_variance=FALSE, max_iter=5000) {
  stopifnot(is.matrix(R))
  stopifnot(nrow(R) == ncol(R))
  stopifnot(nrow(R) == nrow(df))
  if (is.null(n)) {
    n_vals <- unique(df$N)
    n <- if (length(n_vals) == 1) n_vals[1] else median(df$N, na.rm=TRUE)
  }
  z <- df$BETA / df$SE
  fit <- susie_rss(
    z = z,
    R = R,
    n = n,
    L = L,
    estimate_residual_variance = estimate_residual_variance,
    max_iter = max_iter
  )
  cs_list <- fit$sets$cs
  cs_df <- data.frame()
  if (!is.null(cs_list) && length(cs_list) > 0) {
    for (nm in names(cs_list)) {
      idx <- cs_list[[nm]]
      df_cs <- data.frame(
        CS = nm,
        SNP = df$SNP[idx],
        BP = df$BP[idx],
        P = df$P[idx],
        PIP = fit$pip[idx],
        stringsAsFactors = FALSE
      )
      df_cs <- df_cs[order(df_cs$PIP, decreasing=TRUE), ]
      cs_df <- rbind(cs_df, df_cs)
    }
  }
  list(
    fit = fit,
    cs = cs_df
  )
}

setwd("/Users/c24102394/Desktop/neurobridge/outputs/defined_loci/AD_SCZ/locus_chr11_46732038_47732038")
dir.create("susie_rss", showWarnings=FALSE, recursive=TRUE)
R <- as.matrix(read.table(gzfile("ld_gwas/locus.ld.gz"), header=FALSE))
storage.mode(R) <- "double"

ad <- read.table("gwas_AD.ldorder.tsv",  header=TRUE, sep="\t", stringsAsFactors=FALSE)
scz <- read.table("gwas_SCZ.ldorder.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)

stopifnot(nrow(R) == ncol(R))
stopifnot(nrow(R) == nrow(ad))
stopifnot(nrow(R) == nrow(scz))

res_ad <- run_susie(ad, R, L=1, estimate_residual_variance=FALSE, max_iter=5000)
res_scz <- run_susie(scz, R, L=1, estimate_residual_variance=FALSE, max_iter=5000)

res_ad$cs
res_scz$cs










pip <- res_ad$fit$pip
o <- order(pip, decreasing=TRUE)

cum <- cumsum(pip[o])
k <- which(cum >= 0.95)[1]
idx <- o[1:k]

ad_cs95 <- data.frame(
  SNP = ad$SNP[idx],
  BP  = ad$BP[idx],
  P   = ad$P[idx],
  PIP = pip[idx],
  CUM_PIP = cumsum(pip[idx]),
  stringsAsFactors = FALSE
)

nrow(ad_cs95)
tail(ad_cs95, 1)
head(ad_cs95, 30)





pip <- res_scz$fit$pip
o <- order(pip, decreasing=TRUE)

cum <- cumsum(pip[o])
k <- which(cum >= 0.95)[1]
idx <- o[1:k]

scz_cs95 <- data.frame(
  SNP = scz$SNP[idx],
  BP  = scz$BP[idx],
  P   = scz$P[idx],
  PIP = pip[idx],
  CUM_PIP = cumsum(pip[idx]),
  stringsAsFactors = FALSE
)

nrow(scz_cs95)
tail(scz_cs95, 1)
head(scz_cs95, 30)




length(intersect(ad_cs95$SNP, scz_cs95$SNP))
intersect(ad_cs95$SNP, scz_cs95$SNP)

