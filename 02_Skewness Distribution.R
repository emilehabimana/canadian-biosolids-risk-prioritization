# ============================================================
# PANEL B — ABSOLUTE SKEWNESS DISTRIBUTION (FINAL, FULLY FIXED)
# ============================================================

# ---------- 1. PATHS ----------
data_dir <- "C:/Users/Admin/Documents/Data"
figures_dir <- "C:/Users/Admin/Documents/Data/Figures"

if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)
setwd(data_dir)

# ---------- 2. LIBRARIES ----------
library(ggplot2)
library(moments)
library(showtext)

# ---------- 3. FONT ----------
font_add("Arial", "arial.ttf")
showtext_auto()
showtext_opts(dpi = 300)

# ---------- 4. COLORS ----------
TURQUOISE <- "#40E0D0"
ELECTRIC_CRIMSON <- "#FF3F00"
DARK_GRAY <- "#555555"
BLACK <- "#000000"

# ---------- 5. LOAD DATA (CLEAN HEADERS) ----------
raw_data <- read.csv(
  "ec_concentration.csv",
  stringsAsFactors = FALSE,
  na.strings = c("Not found", "NA", ""),
  check.names = TRUE
)

# Remove empty column names
raw_data <- raw_data[, names(raw_data) != "" & !is.na(names(raw_data))]

compound_col <- names(raw_data)[1]
sample_cols <- names(raw_data)[-1]

# Convert concentrations to numeric
raw_data[sample_cols] <- suppressWarnings(
  lapply(raw_data[sample_cols], as.numeric)
)

# ---------- 6. ABSOLUTE SKEWNESS PER COMPOUND (ROBUST) ----------
skewness_values <- numeric()

for (i in seq_len(nrow(raw_data))) {
  vals <- unlist(raw_data[i, sample_cols], use.names = FALSE)
  vals <- vals[is.finite(vals)]
  if (length(vals) >= 3) {
    skewness_values <- c(skewness_values, abs(skewness(vals)))
  }
}

n_total <- length(skewness_values)
cat("✓ Skewness calculated for", n_total, "compounds\n")

# ---------- 7. HISTOGRAM BINS (INCLUDING ZERO BINS) ----------
bin_width <- 0.25
x_max <- ceiling(max(skewness_values)) + bin_width

bin_edges <- seq(0, x_max, by = bin_width)
bin_centers <- (bin_edges[-1] + bin_edges[-length(bin_edges)]) / 2

binned <- cut(
  skewness_values,
  breaks = bin_edges,
  include.lowest = TRUE
)

bin_counts <- table(binned)

bar_data <- data.frame(
  x = bin_centers,
  count = as.numeric(bin_counts)
)

# ---------- 8. DENSITY (SCALED TO COUNTS) ----------
density_obj <- density(skewness_values, bw = "nrd0")
density_scaled <- data.frame(
  x = density_obj$x,
  y = density_obj$y * n_total * bin_width
)
density_scaled <- density_scaled[
  density_scaled$x >= 0 & density_scaled$x <= x_max, ]

# ---------- 9. SUMMARY STATS ----------
median_skew <- round(median(skewness_values), 2)
mean_skew   <- round(mean(skewness_values), 2)
max_count   <- max(bar_data$count)

# Boxplot placement
boxplot_y <- max_count * 0.25
box_width <- max_count * 0.15

# ---------- 10. BUILD PLOT ----------
p <- ggplot() +
  
  # Histogram
  geom_col(
    data = bar_data,
    aes(x = x, y = count),
    fill = TURQUOISE,
    alpha = 0.6,
    color = BLACK,
    linewidth = 0.6,
    width = bin_width * 0.85
  ) +
  
  # Density curve
  geom_line(
    data = density_scaled,
    aes(x = x, y = y),
    color = ELECTRIC_CRIMSON,
    linewidth = 1.2
  ) +
  
  # Horizontal boxplot (whisker box)
  geom_boxplot(
    data = data.frame(Skewness = skewness_values),
    aes(x = Skewness, y = boxplot_y),
    width = box_width,
    fill = DARK_GRAY,
    alpha = 0.7,
    color = BLACK,
    outlier.shape = NA,
    linewidth = 1.2
  ) +
  
  # Extreme skewness threshold
  geom_vline(
    xintercept = 2,
    linetype = "dashed",
    color = ELECTRIC_CRIMSON,
    linewidth = 1.2
  ) +
  
  # Bin counts (bold)
  geom_text(
    data = bar_data,
    aes(x = x, y = -0.6, label = count),
    size = 6,
    fontface = "bold",
    family = "Arial",
    color = BLACK
  ) +
  
  # Mean & Median annotations (26 pt)
  annotate(
    "text",
    x = median_skew,
    y = max_count * 0.85,
    label = paste("Median =", median_skew),
    size = 26 / .pt,
    fontface = "bold",
    family = "Arial",
    color = BLACK
  ) +
  annotate(
    "text",
    x = mean_skew,
    y = max_count * 0.75,
    label = paste("Mean =", mean_skew),
    size = 26 / .pt,
    fontface = "bold",
    family = "Arial",
    color = BLACK
  ) +
  
  # Axis labels
  labs(
    x = expression(paste("Absolute Skewness (|g"[1], "|)")),
    y = "Count (# analytes)"
  ) +
  
  theme_bw() +
  theme(
    plot.title = element_blank(),
    
    axis.title = element_text(
      size = 30,
      face = "bold",
      family = "Arial",
      color = BLACK
    ),
    
    axis.text = element_text(
      size = 22,
      face = "bold",
      family = "Arial",
      color = BLACK
    ),
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = BLACK, linewidth = 1)
  ) +
  
  scale_x_continuous(
    limits = c(0, x_max),
    breaks = seq(0, ceiling(x_max), by = 1),
    expand = c(0, 0)
  ) +
  ylim(-1.2, max_count * 1.15)

# ---------- 11. SAVE ----------
print(p)

ggsave(
  filename = "PanelB_Skewness_Distribution_Final.png",
  plot = p,
  path = figures_dir,
  width = 14,
  height = 10,
  dpi = 300,
  units = "in"
)

cat("✓ Final skewness distribution saved successfully\n")
