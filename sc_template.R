#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(progress)
  library(data.table)
  library(remotes)
  library(phenoscanner)
  library(TwoSampleMR)
  library(MRPRESSO)
  library(mr.raps)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript run-mr.R exposure.tsv outcome.tsv exposure_label outcome_label p_threshold ldlink_token", call. = FALSE)
}

exp_file <- args[1]
out_file <- args[2]
exp_label <- args[3]
out_label <- args[4]
p_threshold <- as.numeric(args[5])
ld_token <- args[6]

outdir <- file.path("outputs/MR/", paste0(exp_label, "_", out_label))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

read_exposure <- function(file, exposure_name) {
  d <- suppressMessages(
    read_exposure_data(
      filename = file,
      sep = "\t",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      eaf_col = "FRQ",
      pval_col = "P",
      samplesize_col = "N",
      phenotype_col = "exposure",
      min_pval = 1e-200
    )
  )
  d$exposure <- exposure_name
  d
}

read_outcome <- function(file, outcome_name) {
  d <- suppressMessages(
    read_outcome_data(
      filename = file,
      sep = "\t",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col = "A1",
      other_allele_col = "A2",
      eaf_col = "FRQ",
      pval_col = "P",
      samplesize_col = "N",
      phenotype_col = "outcome",
      min_pval = 1e-200
    )
  )
  d$outcome <- outcome_name
  d
}

ld_clump_ldlink <- function(dat, clump_kb = 1000, clump_r2 = 0.001) {
  total <- nrow(dat)
  cat(sprintf("[3/6] LD clumping: %d SNPs before clumping\n", total))
  res <- clump_data(
    dat,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_p1 = 1,
    clump_p2 = 1
  )
  kept <- nrow(res)
  cat(sprintf("[3/6] LD clumping: %d SNPs after clumping\n", kept))
  res
}

filter_phenoscanner_all <- function(dat, exp_name, out_name, ps_p = 5e-8) {
  snps <- unique(dat$SNP)
  if (length(snps) == 0) return(dat)
  chunk_size <- 10L
  n_chunks <- ceiling(length(snps) / chunk_size)
  res_list <- vector("list", n_chunks)
  
  for (i in seq_len(n_chunks)) {
    idx <- ((i - 1L) * chunk_size + 1L):min(i * chunk_size, length(snps))
    this_snps <- snps[idx]
    cat(sprintf("\r[4/6] Phenoscanner all-phenotype scan: batch %d/%d", i, n_chunks))
    flush.console()
    r <- tryCatch(
      phenoscanner(snp = this_snps),
      error = function(e) NULL
    )
    if (!is.null(r) && !is.null(r$results)) {
      res_list[[i]] <- r$results
    }
  }
  
  cat("\n")
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]
  if (length(res_list) == 0) return(dat)
  
  hits <- rbindlist(res_list, fill = TRUE)
  if (!("p" %in% names(hits))) {
    pcol <- intersect(names(hits), c("pvalue", "p_value", "P"))
    if (length(pcol) == 0) return(dat)
    setnames(hits, pcol[1], "p")
  }
  if (!("trait" %in% names(hits))) return(dat)
  if (!("snp" %in% names(hits))) return(dat)
  
  hits <- hits[!is.na(p)]
  hits <- hits[p < ps_p]
  if (nrow(hits) == 0) return(dat)
  
  keep_traits <- tolower(c(exp_name, out_name))
  bad <- unique(hits[!(tolower(trait) %in% keep_traits), snp])
  dat[!(dat$SNP %in% bad), ]
}

run_one_direction <- function(exp_file, out_file, prefix, exp_name, out_name, p_threshold, outdir, token) {
  
  cat("=== Running", exp_name, "->", out_name, "===\n")
  
  cat("[1/6] Loading exposure and p-value filtering\n")
  exp <- read_exposure(exp_file, exp_name)
  exp <- exp[exp$pval.exposure < p_threshold, ]
  if (nrow(exp) == 0) {
    cat("No instruments at p <", p_threshold, "for", exp_name, "\n")
    return(invisible())
  }
  fwrite(as.data.table(exp), file.path(outdir, paste0(prefix, "_instruments_raw.csv")))
  
  cat("[2/6] F-stat calculation and weak instrument filter\n")
  exp$F <- (exp$beta.exposure^2) / (exp$se.exposure^2)
  exp$F_pval <- pchisq(exp$F, df = 1, lower.tail = FALSE)
  exp <- exp[exp$F > 10, ]
  if (nrow(exp) == 0) {
    cat("All instruments for", exp_name, "are weak (F <= 10)\n")
    return(invisible())
  }
  fwrite(
    as.data.table(exp[, c("SNP", "beta.exposure", "se.exposure", "pval.exposure", "F", "F_pval")]),
    file.path(outdir, paste0(prefix, "_instruments_with_F.csv"))
  )
  
  cat("[3/6] LD clumping\n")
  exp <- ld_clump_ldlink(exp)
  if (nrow(exp) == 0) {
    cat("No LD-independent instruments remain for", exp_name, "\n")
    return(invisible())
  }
  fwrite(as.data.table(exp), file.path(outdir, paste0(prefix, "_instruments_LDclumped.csv")))
  
  cat("[4/6] Phenoscanner all-phenotype filtering\n")
  exp <- filter_phenoscanner_all(exp, exp_name, out_name, ps_p = 5e-8)
  if (nrow(exp) == 0) {
    cat("All LD-independent instruments flagged by Phenoscanner for", exp_name, "\n")
    return(invisible())
  }
  fwrite(as.data.table(exp), file.path(outdir, paste0(prefix, "_instruments_no_phenoscanner_hits.csv")))
  
  cat("[5/6] Harmonising exposure and outcome\n")
  out <- read_outcome(out_file, out_name)
  out <- out[out$SNP %in% exp$SNP, ]
  harm <- harmonise_data(exp, out)
  if (nrow(harm) == 0) {
    cat("No harmonised SNPs for", exp_name, "->", out_name, "\n")
    return(invisible())
  }
  fwrite(as.data.table(harm), file.path(outdir, paste0(prefix, "_harmonised.csv")))
  
  cat("[6/6] Running MR models and sensitivity analyses (incl. MR-PRESSO)\n")
  
  mr_res <- mr(harm, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
  mr_res$lo_ci <- mr_res$b - 1.96 * mr_res$se
  mr_res$up_ci <- mr_res$b + 1.96 * mr_res$se
  mr_res$OR <- exp(mr_res$b)
  mr_res$OR_lo <- exp(mr_res$lo_ci)
  mr_res$OR_up <- exp(mr_res$up_ci)
  
  het_res <- mr_heterogeneity(harm)
  pleio_res <- mr_pleiotropy_test(harm)
  loo_res <- mr_leaveoneout(harm)
  loo_plot <- mr_leaveoneout_plot(loo_res)
  
  fwrite(as.data.table(mr_res), file.path(outdir, paste0(prefix, "_mr_results.csv")))
  fwrite(as.data.table(het_res), file.path(outdir, paste0(prefix, "_heterogeneity.csv")))
  fwrite(as.data.table(pleio_res), file.path(outdir, paste0(prefix, "_pleiotropy.csv")))
  fwrite(as.data.table(loo_res), file.path(outdir, paste0(prefix, "_leaveoneout.csv")))
  
  pdf(file.path(outdir, paste0(prefix, "_leaveoneout_plot.pdf")))
  print(loo_plot[[1]])
  dev.off()
  
  if (nrow(harm) >= 4) {
    
    cat("Running MR-PRESSO on harmonised SNPs\n")
    
    presso_res <- tryCatch(
      mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = harm,
        NbDistribution = 1000,
        SignifThreshold = 0.05
      ),
      error = function(e) {
        message("MR-PRESSO failed: ", e$message)
        NULL
      }
    )
    
    if (!is.null(presso_res)) {
      
      if (!is.null(presso_res$`Main MR results`)) {
        main_tab <- as.data.table(presso_res$`Main MR results`, keep.rownames = TRUE)
        setnames(main_tab, "rn", "model", skip_absent = TRUE)
        
        bcol <- intersect(names(main_tab), c("Causal Estimate", "Causal.Estimate", "Estimate", "b", "Beta"))
        scol <- intersect(names(main_tab), c("Sd", "SE", "StdError", "se"))
        
        if (length(bcol) > 0 && length(scol) > 0) {
          b <- as.numeric(main_tab[[bcol[1]]])
          se <- as.numeric(main_tab[[scol[1]]])
          main_tab$lo_ci <- b - 1.96 * se
          main_tab$up_ci <- b + 1.96 * se
          main_tab$OR <- exp(b)
          main_tab$OR_lo <- exp(main_tab$lo_ci)
          main_tab$OR_up <- exp(main_tab$up_ci)
        }
        
        fwrite(main_tab, file.path(outdir, paste0(prefix, "_mrpresso_main_withOR.csv")))
      }
      
      if (!is.null(presso_res$`Global Test`)) {
        global_tab <- as.data.table(presso_res$`Global Test`)
        fwrite(global_tab, file.path(outdir, paste0(prefix, "_mrpresso_global.csv")))
      }
      
      outlier_snps <- NULL
      if (!is.null(presso_res$`Outlier Test`)) {
        outlier_tab <- as.data.table(presso_res$`Outlier Test`, keep.rownames = TRUE)
        setnames(outlier_tab, "rn", "SNP", skip_absent = TRUE)
        fwrite(outlier_tab, file.path(outdir, paste0(prefix, "_mrpresso_outliers.csv")))
        outlier_snps <- outlier_tab$SNP
      }
      
      if (!is.null(presso_res$`Distortion Test`)) {
        distortion_tab <- as.data.table(presso_res$`Distortion Test`)
        fwrite(distortion_tab, file.path(outdir, paste0(prefix, "_mrpresso_distortion.csv")))
      }
      
      if (!is.null(outlier_snps) && length(outlier_snps) > 0) {
        
        cat("Re-running MR after removing MR-PRESSO outliers\n")
        
        outdir_outlier <- file.path(outdir, "OUTLIER.REMOVEDRESULTS")
        dir.create(outdir_outlier, showWarnings = FALSE, recursive = TRUE)
        
        harm_no_out <- harm[!(harm$SNP %in% outlier_snps), ]
        if (nrow(harm_no_out) > 0) {
          
          fwrite(
            as.data.table(harm_no_out),
            file.path(outdir_outlier, paste0(prefix, "_harmonised_outlier_removed.csv"))
          )
          
          mr_res2 <- mr(harm_no_out, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode", "mr_raps"))
          mr_res2$lo_ci <- mr_res2$b - 1.96 * mr_res2$se
          mr_res2$up_ci <- mr_res2$b + 1.96 * mr_res2$se
          mr_res2$OR <- exp(mr_res2$b)
          mr_res2$OR_lo <- exp(mr_res2$lo_ci)
          mr_res2$OR_up <- exp(mr_res2$up_ci)
          
          het_res2 <- mr_heterogeneity(harm_no_out)
          pleio_res2 <- mr_pleiotropy_test(harm_no_out)
          loo_res2 <- mr_leaveoneout(harm_no_out)
          loo_plot2 <- mr_leaveoneout_plot(loo_res2)
          
          fwrite(as.data.table(mr_res2), file.path(outdir_outlier, paste0(prefix, "_mr_results_outlier_removed.csv")))
          fwrite(as.data.table(het_res2), file.path(outdir_outlier, paste0(prefix, "_heterogeneity_outlier_removed.csv")))
          fwrite(as.data.table(pleio_res2), file.path(outdir_outlier, paste0(prefix, "_pleiotropy_outlier_removed.csv")))
          fwrite(as.data.table(loo_res2), file.path(outdir_outlier, paste0(prefix, "_leaveoneout_outlier_removed.csv")))
          
          pdf(file.path(outdir_outlier, paste0(prefix, "_leaveoneout_plot_outlier_removed.pdf")))
          print(loo_plot2[[1]])
          dev.off()
          
        } else {
          cat("All SNPs removed as outliers; no data left for outlier-removed MR\n")
        }
      }
    }
    
  } else {
    cat("Not enough SNPs for MR-PRESSO (need >= 4 instruments)\n")
  }
  
  cat("=== Finished", exp_name, "->", out_name, "===\n")
}

run_one_direction(exp_file, out_file, "forward_expToOut", exp_label, out_label, p_threshold, outdir, ld_token)
run_one_direction(out_file, exp_file, "reverse_outToExp", out_label, exp_label, p_threshold, outdir, ld_token)

cat("Done\n")
