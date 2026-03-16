#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(LAVA)
  library(progress)
  library(data.table)
})

options(datatable.fread.nThread = 1)
options(datatable.fread.mmap = FALSE)
data.table::setDTthreads(1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 15) stop("usage: lava_pair.R ref_prefix loci_file info_tsv overlap_csv out_dir ph1 ph1_cases ph1_ctrl ph1_prev ph1_file ph2 ph2_cases ph2_ctrl ph2_prev ph2_file")

ref_prefix <- args[1]
loci_file  <- args[2]
info_tsv   <- args[3]
ov_file    <- args[4]
out_dir    <- args[5]

ph1_name  <- args[6]
ph1_cases <- as.numeric(args[7])
ph1_ctrls <- as.numeric(args[8])
ph1_prev  <- as.numeric(args[9])
ph1_file  <- args[10]

ph2_name  <- args[11]
ph2_cases <- as.numeric(args[12])
ph2_ctrls <- as.numeric(args[13])
ph2_prev  <- as.numeric(args[14])
ph2_file  <- args[15]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

PHENOS   <- c(ph1_name, ph2_name)
in_files <- c(ph1_file, ph2_file)

info <- data.frame(
  phenotype   = PHENOS,
  cases       = c(ph1_cases, ph2_cases),
  controls    = c(ph1_ctrls, ph2_ctrls),
  prevalence  = c(ph1_prev,  ph2_prev),
  filename    = in_files,
  stringsAsFactors = FALSE
)

dir.create(dirname(info_tsv), recursive = TRUE, showWarnings = FALSE)
fwrite(info, info_tsv, sep = "\t", quote = FALSE, na = "NA")

loci   <- read.loci(loci_file)
p_gate <- 0.05 / nrow(loci)

is_mhc <- function(chr, start, stop) {
  chr6  <- as.character(chr) %in% c("6", "chr6")
  inwin <- !(stop < 25000000 | start > 34000000)
  isTRUE(chr6 & inwin)
}

ov <- as.matrix(read.csv(ov_file, row.names = 1, check.names = FALSE))

OV <- diag(2)
dimnames(OV) <- list(PHENOS, PHENOS)
OV[ph1_name, ph2_name] <- OV[ph2_name, ph1_name] <- ov[ph1_name, ph2_name]

overlap_all_file <- file.path(out_dir, paste0("overlap_corr_for_LAVA_", paste(PHENOS, collapse = "_"), ".csv"))
write.table(OV, overlap_all_file, sep = "\t", quote = FALSE, col.names = NA)

inp <- process.input(
  input.info.file      = info_tsv,
  sample.overlap.file  = overlap_all_file,
  ref.prefix           = ref_prefix,
  phenos               = PHENOS
)

pb <- progress_bar$new(total = nrow(loci), width = 60)
rows_bivar <- vector("list", nrow(loci))
ri <- 0L

for (i in seq_len(nrow(loci))) {
  pb$tick()
  loc <- process.locus(loci[i, ], inp)
  
  ph <- intersect(loc$phenos, PHENOS)
  if (length(ph) < 2) next
  
  uni <- run.univ(loc, phenos = ph)
  
  h2 <- setNames(rep(NA_real_, length(PHENOS)), PHENOS)
  pu <- setNames(rep(NA_real_, length(PHENOS)), PHENOS)
  if (!is.null(uni) && nrow(uni) > 0) {
    for (phn in ph) {
      h2[phn] <- uni$h2.obs[uni$phen == phn]
      pu[phn] <- uni$p[uni$phen == phn]
    }
  }
  
  gate1 <- is.finite(h2[ph1_name]) && h2[ph1_name] > 0 && is.finite(pu[ph1_name]) && pu[ph1_name] < p_gate
  gate2 <- is.finite(h2[ph2_name]) && h2[ph2_name] > 0 && is.finite(pu[ph2_name]) && pu[ph2_name] < p_gate
  if (!(gate1 && gate2)) next
  
  bv <- try(run.bivar(loc, phenos = c(ph1_name, ph2_name), param.lim = 2), silent = TRUE)
  if (inherits(bv, "try-error") || is.null(bv) || nrow(bv) == 0) next
  
  lc <- if (
    is.finite(bv$rho) &&
    is.finite(h2[ph1_name]) && is.finite(h2[ph2_name]) &&
    h2[ph1_name] >= 0 && h2[ph2_name] >= 0
  ) {
    bv$rho * sqrt(h2[ph1_name] * h2[ph2_name])
  } else {
    NA_real_
  }
  
  ri <- ri + 1L
  rows_bivar[[ri]] <- data.frame(
    locus = loc$id, chr = loc$chr, start = loc$start, stop = loc$stop,
    pair  = paste(ph1_name, ph2_name, sep = "_"),
    phen1 = ph1_name, phen2 = ph2_name,
    n_snps = if (!is.null(loc$n.snps)) as.integer(loc$n.snps) else NA_integer_,
    h2_1 = h2[ph1_name], p_1 = pu[ph1_name],
    h2_2 = h2[ph2_name], p_2 = pu[ph2_name],
    rho = bv$rho, rho.lower = bv$rho.lower, rho.upper = bv$rho.upper,
    r2 = bv$r2, p = bv$p,
    local_cov = lc,
    stringsAsFactors = FALSE
  )
}

bivar <- rbindlist(rows_bivar[seq_len(ri)], use.names = TRUE, fill = TRUE)

if (nrow(bivar) > 0) {
  bivar[, mhc := is_mhc(as.integer(chr), start, stop)]
  idx <- which(is.finite(bivar$p))
  K <- length(idx)
  if (K > 0) {
    alpha <- 0.05 / K
    bivar$p_bonf_paper <- NA_real_
    bivar$sig_paper    <- NA
    bivar$q_fdr        <- NA_real_
    bivar$p_bonf_paper[idx] <- bivar$p[idx] * K
    bivar$sig_paper[idx]    <- bivar$p[idx] < alpha
    bivar$q_fdr[idx]        <- p.adjust(bivar$p[idx], method = "fdr")
  }
}

fwrite(bivar, file.path(out_dir, "LAVA_local_rg_bivariate.tsv"), sep = "\t")
cat("Saved:\n")
cat(file.path(out_dir, "LAVA_local_rg_bivariate.tsv"), "\n")
cat(file.path(out_dir, "overlap_corr_for_LAVA_", paste(PHENOS, collapse = "_"), ".csv"), "\n")