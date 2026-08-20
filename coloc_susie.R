#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(arrow)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop("Usage: Rscript coloc_susie.R <pqtl_dataset> <protein_id> <pheno_id> <gwas_parquet> <pqtl_parquet> <n_cases> <n_controls> [results_dir]")
}

# same exposure/outcome convention as coloc.R - exposure pQTL == quant ALWAYS
# LD matrices come from bin/ld_matrix.sh, saved next to pqtl_file/gwas_file as
# {protein_id}_ld_matrix.{ld,snplist} and {pheno_id}_ld_matrix.{ld,snplist}

pqtl_dataset <- args[1]
protein_id   <- args[2]
pheno_id     <- args[3]
gwas_file    <- args[4]
pqtl_file    <- args[5]
n_cases      <- as.numeric(args[6])
n_controls   <- as.numeric(args[7])
results_dir  <- ifelse(length(args) >= 8, args[8], "results")

exposure_def <- "quant"
outcome_def  <- "cc"
pp4_thresh   <- 0.70
out_dir <- file.path(results_dir, "coloc_susie", pqtl_dataset)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, paste0(pheno_id, "_", protein_id, "_coloc_susie.tsv"))
protein  <- as.data.table(read_parquet(pqtl_file))
gwas     <- as.data.table(read_parquet(gwas_file))

protein <- protein[order(P)][!duplicated(SNP)]
gwas    <- gwas[order(P)][!duplicated(SNP)]
top_snp <- protein$SNP[1]

n_protein <- max(protein$N)
n_gwas <- n_cases + n_controls
s_gwas <- n_cases / n_gwas

print(paste0("[TRACKING] protein_id: ", protein_id))
print(paste0("[TRACKING] pQTL SNPs: ", nrow(protein), " | GWAS SNPs: ", nrow(gwas)))

# deCODE's own pQTL extraction has used both "FRQ" and "MAF" for the same
# allele-frequency column - pick whichever is actually present
protein_maf_col <- if ("FRQ" %in% names(protein)) "FRQ" else "MAF"

read_ld <- function(prefix) {
  snps <- readLines(paste0(prefix, ".snplist"))
  ld <- as.matrix(fread(paste0(prefix, ".ld"), header = FALSE))
  dimnames(ld) <- list(snps, snps)
  ld
}

ld_dir <- dirname(pqtl_file)
pqtl_ld <- read_ld(file.path(ld_dir, paste0(protein_id, "_ld_matrix")))
gwas_ld <- read_ld(file.path(ld_dir, paste0(pheno_id, "_ld_matrix")))

# restrict + reorder each dataset to the SNPs plink kept in its own LD matrix
# (--maf inside ld_matrix.sh can drop some) - runsusie() needs dataset$snp to
# line up 1:1 with dataset$LD's dimnames
protein <- protein[match(rownames(pqtl_ld), SNP)]
gwas    <- gwas[match(rownames(gwas_ld), SNP)]

dataset1 <- list(
  snp     = protein$SNP,
  beta    = protein$BETA,
  varbeta = protein$SE^2,
  MAF     = protein[[protein_maf_col]],
  N       = n_protein,
  type    = exposure_def,
  LD      = pqtl_ld
)

dataset2 <- list(
  snp     = gwas$SNP,
  beta    = gwas$BETA,
  varbeta = gwas$SE^2,
  MAF     = gwas$FRQ,
  N       = n_gwas,
  s       = s_gwas,
  type    = outcome_def,
  LD      = gwas_ld
)

# SuSiE fine-mapping per trait, then pairwise coloc across credible sets -
# resolves allelic heterogeneity within a cis-region, unlike coloc.abf's
# single causal variant assumption
susie_pqtl <- runsusie(dataset1, suffix = 1)
susie_gwas <- runsusie(dataset2, suffix = 2)

n_pqtl_cs <- length(susie_pqtl$sets$cs)
n_gwas_cs <- length(susie_gwas$sets$cs)
print(paste0("[TRACKING] credible sets - pQTL: ", n_pqtl_cs, " | GWAS: ", n_gwas_cs))

if (n_pqtl_cs == 0 || n_gwas_cs == 0) {
  print(paste0("[CONCERN] No credible sets for ", protein_id, " - writing null result"))
  summary <- data.table(
    nsnps = NA_integer_, hit1 = NA_character_, hit2 = NA_character_,
    PP.H0.abf = NA_real_, PP.H1.abf = NA_real_, PP.H2.abf = NA_real_,
    PP.H3.abf = NA_real_, PP.H4.abf = NA_real_
  )
} else {
  res <- coloc.susie(susie_pqtl, susie_gwas, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  print(res$summary)
  summary <- as.data.table(res$summary)
}

setnames(summary, c("hit1", "hit2"), c("pqtl_hit_snp", "gwas_hit_snp"), skip_absent = TRUE)
summary[, protein_id := protein_id]
summary[, outcome_trait := pheno_id]
summary[, top_snp := top_snp]
summary[, n_pqtl_credible_sets := n_pqtl_cs]
summary[, n_gwas_credible_sets := n_gwas_cs]
summary[, n_pqtl_snps := nrow(protein)]
summary[, n_gwas_snps := nrow(gwas)]
summary[, n_cases := n_cases]
summary[, n_controls := n_controls]
summary[, n_gwas := n_gwas]
summary[, s_gwas := s_gwas]
summary[, pp4_threshold := pp4_thresh]
summary[, coloc_pass := !is.na(PP.H4.abf) & PP.H4.abf >= pp4_thresh]

setcolorder(summary, c(
  "protein_id",
  "outcome_trait",
  "top_snp",
  "pqtl_hit_snp",
  "gwas_hit_snp",
  "nsnps",
  "PP.H0.abf",
  "PP.H1.abf",
  "PP.H2.abf",
  "PP.H3.abf",
  "PP.H4.abf",
  "coloc_pass",
  "n_pqtl_credible_sets",
  "n_gwas_credible_sets",
  "n_pqtl_snps",
  "n_gwas_snps",
  "n_cases",
  "n_controls",
  "n_gwas",
  "s_gwas",
  "pp4_threshold"
))

fwrite(summary, out_file, sep = "\t")
print(paste0("[DONE] Saved COLOC.SUSIE result: ", out_file))
