#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(coloc) })

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 7) stop("Usage: Rscript coloc.R <file1> <file2> <out_file> <type1> <type2> <label1> <label2> [s1] [s2]")

file1 <- args[1]
file2 <- args[2]
out_file <- args[3]
type1 <- args[4]
type2 <- args[5]
label1 <- args[6]
label2 <- args[7]
s1 <- if (length(args) >= 8) as.numeric(args[8]) else 0.5
s2 <- if (length(args) >= 9) as.numeric(args[9]) else 0.5

getN <- function(mm, suff) {
  cand <- c(paste0("N",suff), paste0("N_eff",suff), paste0("Neff",suff), paste0("NEFF",suff))
  cand <- unique(cand)
  for (nm in cand) if (nm %in% names(mm)) return(mm[[nm]])
  rep(NA_real_, nrow(mm))
}

d1 <- read.table(file1, sep="\t", header=TRUE, stringsAsFactors=FALSE)
d2 <- read.table(file2, sep="\t", header=TRUE, stringsAsFactors=FALSE)
m <- merge(d1, d2, by="SNP", suffixes=c(".1",".2"))

if (nrow(m) < 10) {
  out <- data.frame(
    file1=basename(file1),
    file2=basename(file2),
    label1=label1,
    label2=label2,
    PP.H0.abf=NA, PP.H1.abf=NA, PP.H2.abf=NA, PP.H3.abf=NA, PP.H4.abf=NA,
    lead_snp=NA,
    n_snps=nrow(m)
  )
  write.table(out, out_file, sep="\t", quote=FALSE, row.names=FALSE)
  quit(save="no")
}

same <- m$A1.1 == m$A1.2 & m$A2.1 == m$A2.2
flip <- m$A1.1 == m$A2.2 & m$A2.1 == m$A1.2
m$BETA.2[flip] <- -m$BETA.2[flip]
m <- m[same | flip, ]

if (nrow(m) < 10) {
  out <- data.frame(
    file1=basename(file1),
    file2=basename(file2),
    label1=label1,
    label2=label2,
    PP.H0.abf=NA, PP.H1.abf=NA, PP.H2.abf=NA, PP.H3.abf=NA, PP.H4.abf=NA,
    lead_snp=NA,
    n_snps=nrow(m)
  )
  write.table(out, out_file, sep="\t", quote=FALSE, row.names=FALSE)
  quit(save="no")
}

lead_snp <- m$SNP[which.min(m$P.1)]
N1 <- unique(na.omit(getN(m, ".1")))
N2 <- unique(na.omit(getN(m, ".2")))
N1 <- if (length(N1)) as.numeric(N1[1]) else NA_real_
N2 <- if (length(N2)) as.numeric(N2[1]) else NA_real_
ds1 <- list(snp=m$SNP, beta=m$BETA.1, varbeta=m$SE.1^2, N=N1, type=type1)
ds2 <- list(snp=m$SNP, beta=m$BETA.2, varbeta=m$SE.2^2, N=N2, type=type2)

if (type1 == "cc") ds1$s <- s1
if (type2 == "cc") ds2$s <- s2
if (type1 == "quant") ds1$sdY <- 1
if (type2 == "quant") ds2$sdY <- 1

co <- coloc.abf(ds1, ds2)
out <- as.data.frame(t(co$summary), stringsAsFactors=FALSE)
out$file1 <- basename(file1)
out$file2 <- basename(file2)
out$label1 <- label1
out$label2 <- label2
out$lead_snp <- lead_snp
out$n_snps <- nrow(m)

write.table(out, out_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("Wrote:", out_file, "\n")