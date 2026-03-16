#!/usr/bin/env Rscript
# Modular script for running global HDL-L (pairwise)

suppressPackageStartupMessages({
  library(data.table)
  library(HDL)
})

# remotes::install_github("zhenin/HDL", subdir="HDL", upgrade="never")

packageDescription("HDL")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 12) {
  stop(paste(
    "Usage:",
    "hdll_global.R <gwas1_ldsc.tsv> <gwas2_ldsc.tsv> <trait1> <trait2> <beta1_col> <se1_col> <beta2_col> <se2_col> <LD.path> <bim.path> <out_dir> <cov_min>",
    sep = "\n"
  ))
}

gwas1_in <- args[1]
gwas2_in <- args[2]
trait1 <- args[3]
trait2 <- args[4]
beta1_col <- args[5]
se1_col <- args[6]
beta2_col <- args[7]
se2_col <- args[8]
LD.path <- args[9]
bim.path <- args[10]
out_dir <- args[11]
cov_min <- as.numeric(args[12])

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(gwas1_in))
stopifnot(file.exists(gwas2_in))
stopifnot(file.exists(file.path(LD.path, "HDLL_LOC_snps.RData")))

load(file.path(LD.path, "HDLL_LOC_snps.RData"))

compute_z <- function(df, beta_col, se_col) {
  df <- as.data.table(df)
  df[, Z := get(beta_col) / get(se_col)]
  df
}

to_hdll <- function(df, trait, beta_col, se_col) {
  df <- as.data.table(df)
  need <- c("SNP", "A1", "A2", beta_col, se_col, "N")
  miss <- setdiff(need, names(df))
  if (length(miss) > 0) stop(paste(trait, "missing columns:", paste(miss, collapse = ", ")))
  
  df <- df[!is.na(SNP) & !is.na(A1) & !is.na(A2)]
  df <- df[!is.na(get(beta_col)) & !is.na(get(se_col)) & !is.na(N)]
  
  df[, N := as.numeric(N)]
  df[, (beta_col) := as.numeric(get(beta_col))]
  df[, (se_col) := as.numeric(get(se_col))]
  
  df <- df[is.finite(N) & is.finite(get(beta_col)) & is.finite(get(se_col)) & get(se_col) > 0]
  
  df <- compute_z(df, beta_col, se_col)
  df <- df[is.finite(Z)]
  
  df[, SNP := as.character(SNP)]
  df[, A1 := toupper(as.character(A1))]
  df[, A2 := toupper(as.character(A2))]
  
  df <- df[!duplicated(SNP)]
  df <- df[, .(SNP, A1, A2, N, Z)]
  df
}

gwas1_raw <- fread(gwas1_in, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
gwas2_raw <- fread(gwas2_in, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)

gwas1 <- to_hdll(gwas1_raw, trait1, beta1_col, se1_col)
gwas2 <- to_hdll(gwas2_raw, trait2, beta2_col, se2_col)

fwrite(gwas1, file.path(out_dir, paste0(trait1, ".hdll.tsv")), sep = "\t")
fwrite(gwas2, file.path(out_dir, paste0(trait2, ".hdll.tsv")), sep = "\t")

NEWLOC <- as.data.table(NEWLOC)
NEWLOC[, CHR := as.integer(CHR)]
NEWLOC[, piece := as.integer(piece)]

locus_coverage <- function(chr, piece) {
  bim_file <- file.path(bim.path, sprintf("ukb_chr%d.%d_n336000.imputed_clean.bim", chr, piece))
  if (!file.exists(bim_file)) {
    return(data.table(CHR = chr, piece = piece, n_ref = NA_integer_, cov_1 = NA_real_, cov_2 = NA_real_))
  }
  b <- fread(bim_file, header = FALSE, showProgress = FALSE)
  if (ncol(b) < 2) {
    return(data.table(CHR = chr, piece = piece, n_ref = NA_integer_, cov_1 = NA_real_, cov_2 = NA_real_))
  }
  if (ncol(b) >= 6) {
    setnames(b, c("CHR","SNP","CM","BP","A1","A2"))
  } else {
    setnames(b, c("CHR","SNP","CM","BP","A1","A2")[seq_len(ncol(b))])
  }
  b[, SNP := as.character(SNP)]
  n_ref <- nrow(b)
  cov_1 <- mean(b$SNP %chin% gwas1$SNP)
  cov_2 <- mean(b$SNP %chin% gwas2$SNP)
  data.table(CHR = chr, piece = piece, n_ref = n_ref, cov_1 = cov_1, cov_2 = cov_2)
}

cov_list <- lapply(seq_len(nrow(NEWLOC)), function(i) {
  locus_coverage(NEWLOC$CHR[i], NEWLOC$piece[i])
})
cov_tab <- rbindlist(cov_list, fill = TRUE)

fwrite(cov_tab, file.path(out_dir, paste0(trait1, "_", trait2, ".coverage.tsv")), sep = "\t")

idx_all <- seq_len(nrow(NEWLOC))
idx_keep <- idx_all[
  is.finite(cov_tab$cov_1) & is.finite(cov_tab$cov_2) &
    cov_tab$cov_1 >= cov_min & cov_tab$cov_2 >= cov_min
]
if (length(idx_keep) == 0L) idx_keep <- idx_all

runner <- function(chr, piece) {
  HDL.L(
    gwas1 = gwas1, gwas2 = gwas2,
    Trait1name = trait1, Trait2name = trait2,
    LD.path = LD.path, bim.path = bim.path,
    chr = chr, piece = piece,
    eigen.cut = 0.99, N0 = 0, lim = exp(-18), output.file = ""
  )
}

res_list <- lapply(idx_keep, function(i) {
  try(runner(NEWLOC$CHR[i], NEWLOC$piece[i]), silent = TRUE)
})

ok <- vapply(res_list, function(x) !(inherits(x, "try-error") || is.null(x)), logical(1))
if (!any(ok)) {
  fwrite(cov_tab, file.path(out_dir, paste0(trait1, "_", trait2, ".coverage_only.tsv")), sep = "\t")
  stop("No loci returned valid HDL-L results. Check SNP IDs / allele alignment / overlap with bim reference.")
}

res_local <- rbindlist(res_list[ok], fill = TRUE)
if (all(c("chr","piece") %in% names(res_local))) {
  setnames(res_local, c("chr","piece"), c("CHR","piece"))
}

res_join <- merge(res_local, cov_tab, by = c("CHR","piece"), all.x = TRUE)
if (all(c("CHR","piece") %in% names(res_join))) {
  setnames(res_join, c("CHR","piece"), c("chr","piece"))
}

fwrite(res_join, file.path(out_dir, paste0(trait1, "_", trait2, ".local_rg.with_coverage.tsv")), sep = "\t")
saveRDS(res_join, file.path(out_dir, paste0(trait1, "_", trait2, ".local_rg.with_coverage.rds")))

valid <- res_join[
  is.finite(Heritability_1) & is.finite(Heritability_2) & is.finite(Genetic_Covariance) &
    Heritability_1 > 0 & Heritability_2 > 0 &
    is.finite(cov_1) & is.finite(cov_2) &
    cov_1 >= cov_min & cov_2 >= cov_min
]

if (nrow(valid) <= 1L) {
  fwrite(data.table(message = "Not enough valid loci after filtering"),
         file.path(out_dir, paste0(trait1, "_", trait2, ".global_rg.tsv")),
         sep = "\t")
  quit(status = 0)
}

h1 <- valid$Heritability_1
h2 <- valid$Heritability_2
gcov <- valid$Genetic_Covariance

rg_hat <- sum(gcov) / sqrt(sum(h1) * sum(h2))

theta_i <- sapply(seq_len(nrow(valid)), function(i) {
  h1_i <- sum(h1) - h1[i]
  h2_i <- sum(h2) - h2[i]
  gcov_i <- sum(gcov) - gcov[i]
  if (h1_i <= 0 || h2_i <= 0) return(NA_real_)
  gcov_i / sqrt(h1_i * h2_i)
})
theta_i <- theta_i[is.finite(theta_i)]

n <- length(theta_i)
theta_dot <- mean(theta_i)
se_jk <- sqrt(((n - 1) / n) * sum((theta_i - theta_dot)^2))

z <- rg_hat / se_jk
p <- 2 * pnorm(-abs(z))
ci_low <- max(-1, rg_hat - 1.96 * se_jk)
ci_high <- min(1, rg_hat + 1.96 * se_jk)

global_tab <- data.table(
  rg = rg_hat, SE = se_jk, z = z, p = p,
  CI_low = ci_low, CI_high = ci_high,
  n_loci_used = nrow(valid),
  cov_min_used = cov_min
)

fwrite(global_tab, file.path(out_dir, paste0(trait1, "_", trait2, ".global_rg.tsv")), sep = "\t")
saveRDS(global_tab, file.path(out_dir, paste0(trait1, "_", trait2, ".global_rg.rds")))
print(global_tab)