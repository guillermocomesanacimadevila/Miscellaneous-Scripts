library(ggplot2)
library(ggfx)
library(cowplot)

set.seed(1)
df <- data.frame(
  x = rnorm(50),
  y = rnorm(50)
)

p <- ggplot(df, aes(x, y)) +
  geom_point(size = 2.2, alpha = 0.85) +
  labs(title = "Toy scatter", subtitle = "Screenshot-style card") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(18, 18, 18, 18),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

card <- ggdraw() +
  draw_plot(p, x = 0.04, y = 0.04, width = 0.92, height = 0.92) +
  theme(
    plot.background = element_rect(
      fill = "white",
      colour = "#E6E6E6",
      linewidth = 0.6
    )
  )

card_shadow <- with_shadow(
  card,
  sigma = 18,
  x_offset = 8,
  y_offset = -8,
  colour = "black"
)

card_shadow
