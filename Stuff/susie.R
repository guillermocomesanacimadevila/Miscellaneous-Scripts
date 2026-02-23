#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(susieR); library(Matrix) })

setwd("/Users/c24102394/Desktop/neurobridge/outputs/defined_loci/AD_SCZ/locus_chr15_58534174_59534174")

dir.create("susie_rss", showWarnings=FALSE, recursive=TRUE)

R <- as.matrix(read.table(gzfile("ld_gwas/locus.ld.gz"), header=FALSE))
storage.mode(R) <- "double"

ad <- read.table("gwas_AD.ldorder.tsv",  header=TRUE, sep="\t", stringsAsFactors=FALSE)
scz <- read.table("gwas_SCZ.ldorder.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)

stopifnot(nrow(R) == ncol(R))
stopifnot(nrow(R) == nrow(ad))
stopifnot(nrow(R) == nrow(scz))

z_ad <- ad$BETA/ad$SE

colnames(ad)
# extract stuff
bhat <- ad$BETA
shat <- ad$SE
n <- ad$N[1]

fitted_rss <- susie_rss(
  z = z_ad,
  R = R,
  n = n,
  L = 10
)

# summary(fitted_rss)

pip <- fitted_rss$pip
sets <- fitted_rss$sets
cs_list <- sets$cs
stopifnot(!is.null(cs_list), length(cs_list) > 0)

cs_list
for (nm in names(cs_list)) {
  idx <- cs_list[[nm]]
  df_cs <- data.frame(
    CS = nm,
    SNP = ad$SNP[idx],
    BP = ad$BP[idx],
    P = ad$P[idx],
    PIP = pip[idx],
    stringsAsFactors = FALSE
  )
  df_cs <- df_cs[order(df_cs$PIP, decreasing=TRUE), ]
  print(df_cs, row.names = FALSE)
  cat("\n")
}

# check LD
lead_snps <- sapply(fitted_rss$sets$cs, function(idx) ad$SNP[idx])
lead_snps




lead_idx
R_sub <- R[lead_idx, lead_idx]
r2_sub <- R_sub^2
print(lead_snps)
print(round(r2_sub, 3))



fit_ad_L1 <- susie_rss(
  z = z_ad,
  R = R,
  n = n,
  L = 1,
  estimate_residual_variance = FALSE,
  max_iter = 5000
)
fit_ad_L1$sets


cs_idx <- fit_ad_L1$sets$cs[[1]]

ad_cs <- data.frame(
  SNP = ad$SNP[cs_idx],
  BP  = ad$BP[cs_idx],
  P   = ad$P[cs_idx],
  PIP = fit_ad_L1$pip[cs_idx],
  stringsAsFactors = FALSE
)

ad_cs <- ad_cs[order(ad_cs$P), ]
ad_cs

