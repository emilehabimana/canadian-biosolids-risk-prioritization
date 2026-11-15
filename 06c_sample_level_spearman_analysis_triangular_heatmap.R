# ---------- TRIANGULAR + HORIZONTALLY FLIPPED HEATMAP ----------

# Keep only upper triangle and remove diagonal
cor_tri <- cor_mat
cor_tri[lower.tri(cor_tri, diag = TRUE)] <- NA

# Convert to long format
cor_long <- reshape2::melt(
  cor_tri,
  varnames = c("Sample_1", "Sample_2"),
  value.name = "Spearman_rho",
  na.rm = TRUE
)

# Reverse x-axis order for horizontal flip
cor_long$Sample_1 <- factor(
  cor_long$Sample_1,
  levels = rev(unique(cor_long$Sample_1))
)

p <- ggplot(cor_long, aes(Sample_1, Sample_2, fill = Spearman_rho)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "turquoise",
    mid = "red",
    high = "darkblue",
    midpoint = 0,
    limits = c(-1, 1),
    name = expression(rho)
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks  = element_blank(),
    panel.grid  = element_blank()
  ) +
  labs(
    title = "Sample-Level Spearman Correlation of EC Profiles",
    subtitle = "Triangular (unique pairs), horizontally flipped"
  )

ggsave(
  filename = "SampleLevel_Spearman_Correlation_Heatmap_Triangular.png",
  plot = p,
  path = figures_dir,
  width = 14,
  height = 12,
  dpi = 300
)