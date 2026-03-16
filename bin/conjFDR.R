#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(cfdr.pleio)
  library(progress)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: conjFDR.R <input_file> <out_prefix> <ref_dir> <out_dir>")

infile <- args[1]
out_prefix <- args[2]
REF_DIR <- args[3]
OUT_DIR <- args[4]
LOCAL_REF_DIR <- file.path(OUT_DIR, paste0(out_prefix, "_localref"))
PLINK <- Sys.getenv("PLINK", "plink")
REF_BFILE <- Sys.getenv("REF_BFILE", "")
if (REF_BFILE == "") stop("Set REF_BFILE to a merged EUR 1KG reference (PLINK bed/bim/fam prefix)")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Reading input file:", infile, "\n")
pb <- progress_bar$new(
  format = "[:bar] :percent (:current/:total) :elapsed | :message",
  total = 6, clear = FALSE, width = 70
)

pb$tick(0, tokens = list(message = "Loading data..."))

dat <- fread(infile)
pb$tick(tokens = list(message = "Preparing traits..."))
dat[, CHR := as.integer(CHR)]
dat <- dat[CHR >= 1 & CHR <= 22]
dat[, CHR := factor(CHR, levels = 1:22)]

traits <- strsplit(out_prefix, "_", fixed = TRUE)[[1]]
if (length(traits) != 2) stop("out_prefix must look like TRAIT1_TRAIT2, e.g. SCZ_LON")

t1 <- traits[1]
t2 <- traits[2]
p1 <- paste0("P_", t1)
p2 <- paste0("P_", t2)
b1 <- paste0("BETA_", t1)
b2 <- paste0("BETA_", t2)

if (!p1 %in% names(dat)) stop("Missing column: ", p1)
if (!p2 %in% names(dat)) stop("Missing column: ", p2)
if (!b1 %in% names(dat)) stop("Missing column: ", b1)
if (!b2 %in% names(dat)) stop("Missing column: ", b2)

trait1 <- dat[, .(SNP, CHR, POS, BETA = get(b1), PVAL = get(p1))]
trait2 <- dat[, .(SNP, CHR, POS, BETA = get(b2), PVAL = get(p2))]

# trait_scz <- dat[, .(SNP, BETA_SCZ = BETA_SCZ, P_SCZ = P_SCZ)]
# trait_ad  <- dat[, .(SNP, BETA_AD  = BETA_AD,  P_AD  = P_AD)]

pb$tick(tokens = list(message = "Initializing cfdr.pleio..."))
obj <- cfdr_pleio$new()

obj$init_data(
  trait1 = trait1,
  trait2 = trait2,
  trait_names = c(t1, t2),
  refdat = refdata_location(REF_DIR),
  local_refdat_path = LOCAL_REF_DIR,
  verbose = FALSE
)

obj$initialize_pruning_index(n_iter = 50, seed = 154226, verbose = FALSE)
pb$tick(tokens = list(message = "Running conditional FDR computations..."))

obj$calculate_cond_fdr(fdr_trait = 1, verbose = FALSE)
obj$calculate_cond_fdr(fdr_trait = 2, verbose = FALSE)

pb$tick(tokens = list(message = "Extracting and saving results..."))
res <- as.data.table(obj$get_trait_results())

res <- merge(res, dat[, .(SNP, BETA_T1 = get(b1), P_T1 = get(p1))], by="SNP", all.x=TRUE)
res <- merge(res, dat[, .(SNP, BETA_T2 = get(b2), P_T2 = get(p2))], by="SNP", all.x=TRUE)
res[, direction := ifelse(BETA_T1 * BETA_T2 > 0, "Concordant",
                          ifelse(BETA_T1 * BETA_T2 < 0, "Discordant", "Zero"))]

fwrite(res, file.path(OUT_DIR, paste0(out_prefix, "_cfdr_results.tsv")), sep = "\t")
fwrite(res[conj_fdr < 0.05], file.path(OUT_DIR, paste0(out_prefix, "_shared_hits.tsv")), sep = "\t")

fwrite(res[cfdr12 < 0.01], file.path(OUT_DIR, paste0(out_prefix, "_condFDR_", t1, "_given_", t2, "_0.01.tsv")), sep = "\t")
fwrite(res[cfdr21 < 0.01], file.path(OUT_DIR, paste0(out_prefix, "_condFDR_", t2, "_given_", t1, "_0.01.tsv")), sep = "\t")

shared <- res[conj_fdr < 0.05]
if (nrow(shared) > 0) {
  pb$tick(tokens = list(message = "Running PLINK clumping..."))
  
  clump_in <- shared[, .(SNP, P = conj_fdr)]
  clump_file <- file.path(OUT_DIR, paste0(out_prefix, "_clump_input.tsv"))
  fwrite(clump_in, clump_file, sep = "\t")
  
  clump_out_prefix <- file.path(OUT_DIR, paste0(out_prefix, "_clump"))
  cmd <- sprintf('%s --bfile %s --clump %s --clump-field P --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 250 --out %s',
                 PLINK, REF_BFILE, clump_file, clump_out_prefix)
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  clumped_path <- paste0(clump_out_prefix, ".clumped")
  if (file.exists(clumped_path)) {
    cl <- fread(clumped_path, fill = TRUE)
    cl <- cl[!is.na(SNP)]
    lead_snps <- unique(cl$SNP)
    shared_leads <- shared[SNP %in% lead_snps]
    setorder(shared_leads, conj_fdr)
    
    fwrite(shared_leads, file.path(OUT_DIR, paste0(out_prefix, "_shared_leads.tsv")), sep = "\t")
    fwrite(shared_leads[direction == "Concordant"],
           file.path(OUT_DIR, paste0(out_prefix, "_shared_leads_concordant.tsv")), sep = "\t")
    fwrite(shared_leads[direction == "Discordant"],
           file.path(OUT_DIR, paste0(out_prefix, "_shared_leads_discordant.tsv")), sep = "\t")
  }
}

pb$tick(tokens = list(message = "All steps completed successfully."))
cat("\nAll done.\n")
