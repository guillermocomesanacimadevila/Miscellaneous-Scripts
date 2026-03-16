#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(coloc) })

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 7) stop("Usage: Rscript coloc.R <prefix> <trait1> <trait2> <loci_dir> <out_dir> <type1> <type2> [s1] [s2]")

prefix  <- args[1]
trait1  <- args[2]
trait2  <- args[3]
loci_dir<- args[4]
out_dir <- args[5]
type1   <- args[6]
type2   <- args[7]
s1 <- if (length(args) >= 8) as.numeric(args[8]) else 0.5
s2 <- if (length(args) >= 9) as.numeric(args[9]) else 0.5

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
out_file <- file.path(out_dir, paste0(prefix, "_coloc.tsv"))

locus_dirs <- list.dirs(loci_dir, full.names=TRUE, recursive=FALSE)
locus_dirs <- locus_dirs[basename(locus_dirs) != basename(loci_dir)]
if (length(locus_dirs) == 0) locus_dirs <- loci_dir

pairs <- data.frame(
  locus_id = basename(locus_dirs),
  f1 = file.path(locus_dirs, paste0("gwas_", trait1, ".ldgwas.tsv")),
  f2 = file.path(locus_dirs, paste0("gwas_", trait2, ".ldgwas.tsv")),
  stringsAsFactors = FALSE
)

pairs <- pairs[file.exists(pairs$f1) & file.exists(pairs$f2), ]
if (nrow(pairs)==0) stop("No matching locus directories with both trait files between ", trait1, " and ", trait2, ".")

getN <- function(mm, suff) {
  cand <- c(paste0("N",suff), paste0("N_eff",suff), paste0("Neff",suff), paste0("NEFF",suff))
  cand <- unique(cand)
  for (nm in cand) if (nm %in% names(mm)) return(mm[[nm]])
  rep(NA_real_, nrow(mm))
}

results <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {
  locus_id <- pairs$locus_id[i]
  d1 <- read.table(pairs$f1[i], sep="\t", header=TRUE, stringsAsFactors=FALSE)
  d2 <- read.table(pairs$f2[i], sep="\t", header=TRUE, stringsAsFactors=FALSE)
  
  m <- merge(d1, d2, by="SNP", suffixes=c(".1",".2"))
  if (nrow(m) < 10) {
    results[[i]] <- data.frame(PP.H0.abf=NA,PP.H1.abf=NA,PP.H2.abf=NA,PP.H3.abf=NA,PP.H4.abf=NA,nsnps=NA,locus_id=locus_id,lead_snp=NA,n_snps=nrow(m))
    next
  }
  
  same <- m$A1.1==m$A1.2 & m$A2.1==m$A2.2
  flip <- m$A1.1==m$A2.2 & m$A2.1==m$A1.2
  m$BETA.2[flip] <- -m$BETA.2[flip]
  m <- m[same | flip, ]
  
  if (nrow(m) < 10) {
    results[[i]] <- data.frame(PP.H0.abf=NA,PP.H1.abf=NA,PP.H2.abf=NA,PP.H3.abf=NA,PP.H4.abf=NA,nsnps=NA,locus_id=locus_id,lead_snp=NA,n_snps=nrow(m))
    next
  }
  
  lead_snp <- m$SNP[which.min(m$P.1)]
  
  N1 <- unique(na.omit(getN(m,".1")))
  N2 <- unique(na.omit(getN(m,".2")))
  N1 <- if (length(N1)) as.numeric(N1[1]) else NA_real_
  N2 <- if (length(N2)) as.numeric(N2[1]) else NA_real_
  
  ds1 <- list(snp=m$SNP, beta=m$BETA.1, varbeta=m$SE.1^2, N=N1, type=type1)
  ds2 <- list(snp=m$SNP, beta=m$BETA.2, varbeta=m$SE.2^2, N=N2, type=type2)
  
  if (type1=="cc") ds1$s <- s1
  if (type2=="cc") ds2$s <- s2
  if (type1=="quant") ds1$sdY <- 1
  if (type2=="quant") ds2$sdY <- 1
  
  co <- coloc.abf(ds1, ds2)
  s <- as.data.frame(t(co$summary), stringsAsFactors=FALSE)
  s$locus_id <- locus_id
  s$lead_snp <- lead_snp
  s$n_snps <- nrow(m)
  results[[i]] <- s
}

final <- do.call(rbind, results)
final <- final[order(final$locus_id), ]
write.table(final, out_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("Wrote:", out_file, "\n")