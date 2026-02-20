#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(GenomicSEM);library(data.table)})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 14) stop("usage: Rscript ldsc.R pheno1_file pheno2_file pheno1_name pheno2_name case1 ctrl1 case2 ctrl2 pop_prev1 pop_prev2 out_root")

pheno1_file <- args[1]
pheno2_file <- args[2]
pheno1_name <- args[3]
pheno2_name <- args[4]
case1 <- as.numeric(args[5])
ctrl1 <- as.numeric(args[6])
case2 <- as.numeric(args[7])
ctrl2 <- as.numeric(args[8])
pop_prev1 <- as.numeric(args[9])
pop_prev2 <- as.numeric(args[10])
out_root <- args[11]
hm3_path <- args[12]
ld_path <- args[13]
wld_path <- args[14]

#hm3_path <- "/Users/c24102394/ref/ldsc/w_hm3.snplist"
#ld_path  <- "/Users/c24102394/ref/ldsc/eur_w_ld_chr"
#wld_path <- "/Users/c24102394/ref/ldsc/weights_hm3_no_hla"

in_files <- c(pheno1_file, pheno2_file)
trait_names <- c(pheno1_name, pheno2_name)

N1_eff <- 4/(1/case1 + 1/ctrl1)
N2_eff <- 4/(1/case2 + 1/ctrl2)
N_vec <- c(N1_eff, N2_eff)

sample_prev1 <- case1/(case1 + ctrl1)
sample_prev2 <- case2/(case2 + ctrl2)
sample_prev <- c(sample_prev1, sample_prev2)

pop_prev <- c(pop_prev1, pop_prev2)

munge(
  files = in_files,
  maf.filter = 0.01,
  info.filter = 0.9,
  N = N_vec,
  hm3 = hm3_path,
  trait.names = trait_names
)

out_files <- file.path(getwd(), paste0(trait_names, ".sumstats.gz"))
stopifnot(all(file.exists(out_files)))

trait1_qc_dir <- dirname(pheno1_file)
trait1_root <- dirname(trait1_qc_dir)
trait1_ldscdir <- file.path(trait1_root, "post-ldsc")

trait2_qc_dir <- dirname(pheno2_file)
trait2_root <- dirname(trait2_qc_dir)
trait2_ldscdir <- file.path(trait2_root, "post-ldsc")

dir.create(trait1_ldscdir, showWarnings=FALSE, recursive=TRUE)
dir.create(trait2_ldscdir, showWarnings=FALSE, recursive=TRUE)

file.copy(out_files[1], file.path(trait1_ldscdir, paste0(trait_names[1], ".sumstats.gz")), overwrite=TRUE)
file.copy(out_files[2], file.path(trait2_ldscdir, paste0(trait_names[2], ".sumstats.gz")), overwrite=TRUE)

ldsc_out <- ldsc(
  traits = out_files,
  sample.prev = sample_prev,
  population.prev = pop_prev,
  ld = ld_path,
  wld = wld_path,
  ldsc.log = TRUE
)

uni_1 <- ldsc(
  traits = out_files[1],
  ld = ld_path,
  wld = wld_path,
  population.prev = pop_prev[1],
  sample.prev = sample_prev[1]
)

uni_2 <- ldsc(
  traits = out_files[2],
  ld = ld_path,
  wld = wld_path,
  population.prev = pop_prev[2],
  sample.prev = sample_prev[2]
)

get_mean <- function(obj){
  df <- as.data.frame(obj$ldscoutput)
  k <- grep("Mean.*Chi", names(df), ignore.case=TRUE)
  if (length(k)) as.numeric(df[[k[1]]]) else NA_real_
}

qc_out <- rbindlist(list(
  data.table(
    Trait = pheno1_name,
    h2 = as.numeric(uni_1$S[1,1]),
    h2_se = sqrt(as.numeric(uni_1$V[1,1])),
    mean_chisq = get_mean(uni_1),
    Intercept = as.numeric(uni_1$I[1,1])
  ),
  data.table(
    Trait = pheno2_name,
    h2 = as.numeric(uni_2$S[1,1]),
    h2_se = sqrt(as.numeric(uni_2$V[1,1])),
    mean_chisq = get_mean(uni_2),
    Intercept = as.numeric(uni_2$I[1,1])
  )
))

qc_out[, h2_z := h2/h2_se]
qc_out[, pass_h2z := h2_z > 2]
qc_out[, pass_mean := ifelse(is.na(mean_chisq), NA, mean_chisq > 1.02)]
qc_out[, pass_int := Intercept >= 0.9 & Intercept <= 1.1]
qc_out[, pass_all := ifelse(is.na(pass_mean), pass_h2z & pass_int, pass_h2z & pass_mean & pass_int)]
print(qc_out[,.(Trait,h2,h2_se,h2_z,mean_chisq,Intercept,pass_all)])

pair_name <- paste0(pheno1_name,"_",pheno2_name)
outdir <- file.path(out_root, pair_name)
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

fwrite(qc_out, file.path(outdir, paste0("ldsc_trait_qc_",pair_name,".tsv")), sep="\t")
write.csv(ldsc_out$S, file.path(outdir,"ldsc_S.csv"), row.names=FALSE)
write.csv(ldsc_out$V, file.path(outdir,"ldsc_V.csv"), row.names=FALSE)
write.csv(ldsc_out$N, file.path(outdir,"ldsc_N.csv"), row.names=FALSE)
write.csv(data.frame(m=ldsc_out$m), file.path(outdir,"ldsc_m.csv"), row.names=FALSE)

ldsc_rg_summary <- function(ld){
  S <- ld$S
  V <- ld$V
  rg <- S[1,2] / sqrt(S[1,1] * S[2,2])
  g <- c(-0.5*rg/S[1,1], 1/sqrt(S[1,1]*S[2,2]), -0.5*rg/S[2,2])
  se <- sqrt(as.numeric(t(g) %*% V %*% g))
  z <- rg/se
  p <- 2*pnorm(-abs(z))
  zc <- qnorm(0.975)
  ci <- c(max(-1, rg-zc*se), min(1, rg+zc*se))
  data.frame(rg=rg, SE=se, z=z, p=p, CI_low=ci[1], CI_high=ci[2])
}

global <- ldsc_rg_summary(ldsc_out)
print(global)
fwrite(global, file.path(outdir, paste0("ldsc_rg_",pair_name,".tsv")), sep="\t", row.names=FALSE)

I_mat <- ldsc_out$I
rownames(I_mat) <- colnames(I_mat) <- trait_names
overlap_corr <- I_mat / sqrt(outer(diag(I_mat), diag(I_mat), "*"))
diag(overlap_corr) <- 1
write.csv(overlap_corr, file.path(outdir, paste0("overlap_corr_for_LAVA_",pair_name,".csv")), quote=FALSE)
print(overlap_corr)
