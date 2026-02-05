#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(locusplotr)
  library(patchwork)
  library(grid)
})

setwd("/Users/c24102394/Desktop/Desktop - Cardiff University K427JYM6QV/PhD/DiscoveryPipeline/plots/figure3")

lead_snp <- "rs11039131"
sig_y <- -log10(0.05)

locus <- read.table("locus_11_locuszoom.tsv", header = TRUE, sep = "\t")

locus <- locus %>%
  transmute(
    rsid = SNP,
    chromosome = as.integer(CHR),
    position = as.integer(POS),
    effect_allele = A1,
    other_allele = A2,
    p_value = as.numeric(conj_fdr)
  ) %>%
  filter(!is.na(p_value), is.finite(p_value), p_value > 0, p_value <= 1)

p_top <- gg_locusplot(
  df = locus,
  lead_snp = lead_snp,
  rsid = rsid,
  chrom = chromosome,
  pos = position,
  ref = effect_allele,
  alt = other_allele,
  p_value = p_value,
  plot_genes = FALSE,
  plot_recombination = TRUE
) +
  geom_hline(yintercept = sig_y, linetype = "dashed", linewidth = 0.6, colour = "red") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(6, 6, 0, 6),
    panel.spacing = unit(0, "pt")
  )

p_with_genes <- gg_locusplot(
  df = locus,
  lead_snp = lead_snp,
  rsid = rsid,
  chrom = chromosome,
  pos = position,
  ref = effect_allele,
  alt = other_allele,
  p_value = p_value,
  plot_genes = TRUE,
  plot_recombination = FALSE
)

p_genes <- p_with_genes[[2]] +
  theme(
    plot.margin = margin(0, 6, 6, 6),
    panel.spacing = unit(0, "pt")
  )

p_final <- p_top / p_genes +
  plot_layout(heights = c(3.2, 1)) &
  theme(plot.margin = margin(0, 0, 0, 0))

print(p_final)

ggsave("locusplot_conjFDR_sigline_p0.05.png", plot = p_final, width = 10, height = 4, units = "in", dpi = 600)
``
