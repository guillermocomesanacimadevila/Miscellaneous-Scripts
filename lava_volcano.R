#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(grid)
  library(RColorBrewer)
})

centromeres <- c(
  `1`=121535434, `2`=92326171, `3`=90504854, `4`=49660117,
  `5`=46405641, `6`=58830166, `7`=58054331, `8`=43838887,
  `9`=47367679, `10`=39254935, `11`=51644205, `12`=34856694,
  `13`=16000000, `14`=16000000, `15`=17000000, `16`=35335801,
  `17`=22263006, `18`=15460898, `19`=24681782, `20`=26369569,
  `21`=11288129, `22`=13000000
)

results <- read_tsv(
  "/Users/c24102394/Desktop/PhD/AD_SCZ_AGE/outputs/lava/ad_scz_age/LAVA_local_rg_bivariate.tsv",
  show_col_types = FALSE
) %>%
  mutate(
    chr_i = as.integer(chr),
    neglog10_p = -log10(p)
  )

pairs <- unique(results$pair)
shape_vals <- c(16, 15, 17)
markers <- setNames(shape_vals[seq_along(pairs)], pairs)

results <- results %>%
  mutate(pair_lab = gsub("AGE", "LON", pair))

fdr_sig <- results$q_fdr < 0.05
bonf_sig <- as.logical(results$sig_paper)
sig <- fdr_sig | bonf_sig

cap <- 7
results <- results %>%
  mutate(
    above = neglog10_p > cap,
    y_plot = pmin(neglog10_p, cap),
    fdr_sig = fdr_sig,
    bonf_sig = bonf_sig,
    sig = sig
  )

norm01 <- function(x, vmin=-1, vmax=1) pmin(pmax((x - vmin) / (vmax - vmin), 0), 1)
z <- norm01(results$rho, -1, 1)
pal <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
base_cols <- col2rgb(pal(256)[pmin(pmax(floor(z * 255) + 1, 1), 256)], alpha = TRUE) / 255

ycap <- results$y_plot
tval <- 1.08 + 0.22 * sqrt(pmin(pmax(ycap / cap, 0), 1))
tval_mat <- matrix(tval, nrow = 4, ncol = length(tval), byrow = TRUE)
white <- matrix(c(1, 1, 1, 1), nrow = 4, ncol = length(tval), byrow = FALSE)
rgba_mat <- white * (1 - tval_mat) + base_cols * tval_mat
rgba_mat[4, ] <- 1
results$base_rgba <- rgb(rgba_mat[1, ], rgba_mat[2, ], rgba_mat[3, ], rgba_mat[4, ])

dark_blue <- rgb(0.08, 0.12, 0.35, 1.0)
red <- rgb(0.90, 0.32, 0.32, 1.0)
ink <- alpha("black", 0.86)

extreme <- results %>% slice(which.max(neglog10_p))

arm_pos_label <- function(row) {
  mid <- (row$start + row$stop) / 2
  cen <- centromeres[as.character(as.integer(row$chr))]
  arm <- ifelse(mid < cen, "p", "q")
  sprintf("%d%s:%.2f–%.2f", as.integer(row$chr), arm, row$start/1e6, row$stop/1e6)
}

xspan <- max(results$rho, na.rm = TRUE) - min(results$rho, na.rm = TRUE)
dx <- ifelse(is.finite(xspan) && xspan > 0, 0.022 * xspan, 0.03)
dy <- 0.15
dy20 <- 0.36
dy7below <- 0.30
dy4up <- 0.26

sig_to_label <- results %>%
  filter(sig & !above) %>%
  arrange(desc(neglog10_p), p) %>%
  slice_head(n = 25) %>%
  rowwise() %>%
  mutate(
    label = arm_pos_label(cur_data()),
    x = as.numeric(rho),
    y = as.numeric(y_plot),
    is_chr7p_special = (as.integer(chr) == 7) && (start / 1e6 >= 36.49) && (start / 1e6 <= 36.53) && (stop / 1e6 >= 37.96) && (stop / 1e6 <= 38.00),
    is_chr7q_special = (as.integer(chr) == 7) && (abs(start / 1e6 - 93.69) < 0.03) && (abs(stop / 1e6 - 95.17) < 0.03),
    ox = case_when(
      as.integer(chr) == 6 ~ -dx,
      as.integer(chr) == 20 ~ 0.0,
      is_chr7p_special ~ 0.0,
      is_chr7q_special ~ dx * 1.65,
      as.integer(chr) == 4 ~ 0.0,
      as.integer(chr) %in% c(8, 13, 17) ~ dx,
      TRUE ~ ifelse(x >= 0, dx, -dx)
    ),
    oy = case_when(
      as.integer(chr) == 6 ~ dy,
      as.integer(chr) == 20 ~ -dy20,
      is_chr7p_special ~ -dy7below,
      is_chr7q_special ~ dy * 1.38,
      as.integer(chr) == 4 ~ dy4up,
      as.integer(chr) %in% c(8, 13, 17) ~ dy,
      TRUE ~ ifelse(y > (cap - 0.6), -dy, dy)
    ),
    ha = case_when(
      as.integer(chr) == 6 ~ "right",
      as.integer(chr) == 20 ~ "center",
      is_chr7p_special ~ "center",
      is_chr7q_special ~ "left",
      as.integer(chr) == 4 ~ "center",
      as.integer(chr) %in% c(8, 13, 17) ~ "left",
      TRUE ~ ifelse(x >= 0, "left", "right")
    ),
    va = case_when(
      as.integer(chr) == 20 ~ "top",
      is_chr7p_special ~ "top",
      as.integer(chr) == 4 ~ "bottom",
      TRUE ~ "center"
    ),
    label_priority = ifelse(as.integer(chr) == 10, 1L, 0L)
  ) %>%
  ungroup() %>%
  arrange(desc(label_priority), desc(neglog10_p), p)

top_label <- tibble(
  x = as.numeric(extreme$rho),
  y = cap,
  label = arm_pos_label(extreme),
  ox = 0.0,
  oy = 0.28
)

label_box_fill <- alpha("white", 0.86)
label_box_line <- alpha("black", 0.06)

base_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = alpha("black", 0.06), linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title = element_text(color = ink),
    axis.text = element_text(color = alpha("black", 0.72)),
    axis.line = element_line(color = alpha("black", 0.18), linewidth = 0.35),
    legend.position = c(0.985, 0.985),
    legend.justification = c(1, 1),
    legend.title = element_text(face = "bold", size = 10.5, color = ink, margin = margin(b = 4, l = 9.5)),
    legend.text = element_text(size = 9.5, color = alpha("black", 0.80)),
    legend.background = element_rect(fill = alpha("white", 0.90), color = alpha("black", 0.10), linewidth = 0.4),
    legend.box.background = element_rect(fill = alpha("white", 0.90), color = alpha("black", 0.10), linewidth = 0.4),
    legend.key = element_rect(fill = NA, color = NA),
    legend.spacing.y = unit(2.6, "pt"),
    legend.margin = margin(7, 9, 7, 9),
    plot.title = element_blank()
  )

p_main <- ggplot() +
  geom_point(
    data = results %>% filter(!above),
    aes(x = rho, y = y_plot, shape = pair_lab),
    size = 4.5,
    stroke = 0,
    color = results$base_rgba[!results$above],
    alpha = 0.96
  ) +
  geom_point(
    data = results %>% filter(!above, fdr_sig, chr_i != 8),
    aes(x = rho, y = y_plot, shape = pair_lab),
    size = 3.8,
    color = dark_blue,
    stroke = 0.85
  ) +
  geom_point(
    data = results %>% filter(!above, fdr_sig, chr_i == 8),
    aes(x = rho, y = y_plot),
    shape = 8,
    size = 5.0,
    color = dark_blue,
    stroke = 0.9
  ) +
  geom_point(
    data = results %>% filter(!above, bonf_sig, chr_i != 8),
    aes(x = rho, y = y_plot, shape = pair_lab),
    size = 4.9,
    color = dark_blue,
    stroke = 1.05
  ) +
  geom_point(
    data = results %>% filter(!above, bonf_sig, chr_i == 8),
    aes(x = rho, y = y_plot),
    shape = 8,
    size = 6.5,
    color = dark_blue,
    stroke = 1.05
  ) +
  geom_point(
    data = results %>% filter(above),
    aes(x = rho, y = cap, shape = pair_lab),
    size = 5.0,
    color = dark_blue,
    stroke = 1.0
  ) +
  geom_point(
    data = extreme,
    aes(x = rho, y = cap, shape = gsub("AGE", "LON", pair)),
    size = 5.2,
    color = red,
    stroke = 1.15
  ) +
  geom_point(
    data = extreme,
    aes(x = rho, y = cap),
    shape = 8,
    size = 3.0,
    color = "black"
  ) +
  geom_label(
    data = sig_to_label,
    aes(x = x + ox, y = y + oy, label = label),
    size = 2.55,
    fontface = "bold",
    fill = label_box_fill,
    color = ink,
    label.size = 0.25,
    label.r = unit(0.14, "lines"),
    label.padding = unit(0.18, "lines"),
    label.colour = label_box_line,
    hjust = ifelse(sig_to_label$ha == "left", 0, ifelse(sig_to_label$ha == "right", 1, 0.5)),
    vjust = ifelse(sig_to_label$va == "top", 1, ifelse(sig_to_label$va == "bottom", 0, 0.5))
  ) +
  geom_label(
    data = top_label,
    aes(x = x + ox, y = y + oy, label = label),
    size = 2.75,
    fontface = "bold",
    fill = alpha("white", 0.92),
    color = ink,
    label.size = 0.25,
    label.r = unit(0.16, "lines"),
    label.padding = unit(0.20, "lines"),
    label.colour = alpha("black", 0.08),
    hjust = 0.5,
    vjust = 0
  ) +
  scale_shape_manual(
    name = "Trait pair",
    values = setNames(unname(markers), unique(results$pair_lab)),
    guide = guide_legend(
      override.aes = list(
        size = 5.0,
        stroke = 0,
        color = alpha("black", 0.85)
      )
    )
  ) +
  coord_cartesian(ylim = c(0, cap), clip = "off") +
  labs(
    x = expression(Local~r[g]),
    y = expression(-log[10](p))
  ) +
  base_theme

inset_xpad <- 0.06
inset_xlim <- c(as.numeric(extreme$rho) - inset_xpad, as.numeric(extreme$rho) + inset_xpad)
inset_ylim <- c(as.numeric(extreme$neglog10_p) * 0.96, as.numeric(extreme$neglog10_p) * 1.04)

arrow_x0 <- as.numeric(extreme$rho)
arrow_x1 <- inset_xlim[1] + 0.010
arrow_y <- as.numeric(extreme$neglog10_p)

p_inset <- ggplot() +
  geom_point(
    data = extreme,
    aes(x = rho, y = neglog10_p, shape = gsub("AGE", "LON", pair)),
    size = 5.2,
    color = red,
    stroke = 1.1
  ) +
  geom_point(
    data = extreme,
    aes(x = rho, y = neglog10_p),
    shape = 8,
    size = 3.2,
    color = "black"
  ) +
  geom_segment(
    aes(x = arrow_x0, y = arrow_y, xend = arrow_x1, yend = arrow_y),
    linewidth = 0.35,
    color = alpha("black", 0.70),
    arrow = arrow(type = "closed", length = unit(2.4, "mm"))
  ) +
  scale_shape_manual(values = setNames(unname(markers), unique(results$pair_lab))) +
  labs(y = expression(-log[10](p)), x = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = alpha("white", 0.94), color = alpha("black", 0.16), linewidth = 0.6),
    panel.background = element_rect(fill = alpha("white", 0.94), color = NA),
    panel.grid.major = element_line(color = alpha("black", 0.06), linewidth = 0.25),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(color = ink),
    axis.text.y = element_text(color = alpha("black", 0.72)),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(6, 8, 6, 8)
  ) +
  coord_cartesian(xlim = inset_xlim, ylim = inset_ylim)

if (any(results$above)) {
  final_plot <- ggdraw(p_main) +
    draw_plot(p_inset, x = 0.44, y = 0.58, width = 0.30, height = 0.32)
  print(final_plot)
} else {
  print(p_main)
}

out_plot <- if (exists("final_plot")) final_plot else p_main

ggsave(
  filename = "lava_plot.pdf",
  plot = out_plot,
  width = 10.5,
  height = 7,
  units = "in",
  dpi=600
)
