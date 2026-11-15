library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(patchwork)

# -------------------------------
# Paths
# -------------------------------
data_dir    <- "C:/Users/Admin/Documents/Article_2/Data/Data"
figures_dir <- "C:/Users/Admin/Documents/Article_2/Data/Figures"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# Figure size & font settings
# -------------------------------
fig_width_in  <- 90 / 25.4     # 90 mm
fig_height_in <- 1063 / 300    # at least 1063 px @ 300 dpi
text_size     <- 8.5           # SAME size everywhere

# -------------------------------
# EC order (TOP → BOTTOM)
# -------------------------------
ec_order <- c(
  "ROX","CBZ-2OH-2H","TMP","SMZ","ATZ","RNT","ASP","FLO","CBZ-2H",
  "CEL","OMP","FLN","DNT","NFZ","SAL","CLT","CAP","SDZ","GLZ","TAM",
  "MEL-A","STI","AZA","STZ","BZF","ACT","PrP","TCP","TMX"
)

# ============================================================
# HEATMAP – standardized abundance
# ============================================================

ec_ind <- read_csv(
  file.path(data_dir, "ec_robust_scaled_logtrans_MUN_indicators_analysis.csv"),
  show_col_types = FALSE
)

heat_plot_tbl <- ec_ind %>%
  select(Municipality_code, Cluster_NMDS, all_of(ec_order)) %>%
  pivot_longer(
    cols = -c(Municipality_code, Cluster_NMDS),
    names_to = "EC",
    values_to = "StdAbundance"
  ) %>%
  group_by(EC, Cluster_NMDS) %>%
  summarise(
    StdAbundance = mean(StdAbundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    EC = factor(EC, levels = rev(ec_order)),
    Cluster_NMDS = factor(Cluster_NMDS, levels = c(1, 2, 3, 4))
  )

p_heat <- ggplot(
  heat_plot_tbl,
  aes(x = Cluster_NMDS, y = EC, fill = StdAbundance)
) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low  = "#70ffee",
    high = "#fb0101",
    name = "Standardized abundance"
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    axis.text  = element_text(size = text_size, color = "black"),
    axis.title = element_text(size = text_size, color = "black"),
    legend.text  = element_text(size = text_size),
    legend.title = element_text(size = text_size)
  ) +
  labs(
    x = "Municipality cluster",
    y = "Indicator compound"
  )

# ============================================================
# LOLLIPOP – centered, proper legend, same size
# ============================================================

indval_tbl <- read_csv(
  file.path(data_dir, "Indicator_ECs_IndVal_Table.csv"),
  show_col_types = FALSE
)

lollipop_tbl <- indval_tbl %>%
  filter(EC %in% ec_order) %>%
  mutate(
    EC = factor(EC, levels = rev(ec_order)),
    neg_log_p = -log10(p.value)
  )

xmax <- max(c(lollipop_tbl$stat, lollipop_tbl$neg_log_p))

p_lollipop <- ggplot(lollipop_tbl, aes(y = EC)) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_segment(
    aes(x = -neg_log_p, xend = stat),
    color = "grey60"
  ) +
  geom_point(aes(x = stat, color = "IndVal"), size = 3) +
  geom_point(aes(x = -neg_log_p, color = "-log(p-value)"), size = 3) +
  scale_color_manual(
    values = c(
      "IndVal"        = "#fb0101",
      "-log(p-value)" = "#70ffee"
    ),
    name = NULL
  ) +
  coord_cartesian(xlim = c(-xmax, xmax)) +
  scale_x_continuous(
    labels = function(x) ifelse(x == 0, "0", abs(x))
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    axis.text       = element_text(size = text_size, color = "black"),
    axis.title      = element_text(size = text_size, color = "black"),
    legend.text     = element_text(size = text_size),
    plot.title      = element_text(size = text_size),
    legend.position = "top"
  ) +
  labs(
    x = "● IndVal      ● −log(p-value)",
    y = "Indicator compound"
  ) +

# ============================================================
# COMBINE PANELS A / B
# ============================================================

p_final <- p_heat + p_lollipop +
  plot_layout(widths = c(1, 1))

ggsave(
  file.path(figures_dir, "Indicator_Figure_AB.png"),
  plot   = p_final,
  width  = 2 * fig_width_in,   # two panels, each 90 mm
  height = fig_height_in,
  dpi    = 300
)

cat("✓ Figure saved with 90 mm panels and uniform Arial typography\n")