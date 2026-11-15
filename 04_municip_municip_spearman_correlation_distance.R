# ============================================================
# SPEARMAN CORRELATION HEATMAP OF MUNICIPAL EC FINGERPRINTS
# ============================================================

# ---------- PATHS ----------
data_dir    <- "C:/Users/Admin/Documents/Data/Data"
figures_dir <- "C:/Users/Admin/Documents/Data/Figures"
reports_dir <- "C:/Users/Admin/Documents/Data/Reports"

dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(reports_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- LIBRARIES ----------
library(ggplot2)
library(reshape2)
library(rmarkdown)

# ---------- LOAD SPEARMAN CORRELATION DISTANCE MATRIX ----------
dist_df <- read.csv(
  file.path(data_dir, "Municipality_Spearman_Correlation_Distance_MUN.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE,
  row.names = 1
)

# Convert to numeric matrix
dist_mat <- as.matrix(dist_df)
mode(dist_mat) <- "numeric"

# ---------- CONVERT DISTANCE TO CORRELATION ----------
# rho = 1 - distance
cor_mat <- 1 - dist_mat

# ---------- PREPARE DATA FOR HEATMAP ----------
cor_long <- melt(cor_mat, varnames = c("Municipality_1", "Municipality_2"),
                 value.name = "Spearman_rho")

# ---------- PLOT HEATMAP ----------
p <- ggplot(cor_long, aes(Municipality_1, Municipality_2, fill = Spearman_rho)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "turquoise",
    mid = "skyblue",
    high = "darkblue",
    midpoint = 0,
    limits = c(-1, 1),
    name = expression(rho)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Spearman Correlation Heatmap of Municipal EC Fingerprints",
    subtitle = "Monotonic similarity based on robust-scaled EC profiles"
  )

# ---------- SAVE FIGURE ----------
ggsave(
  filename = "Municipality_Spearman_Correlation_Heatmap.png",
  plot = p,
  path = figures_dir,
  width = 12,
  height = 10,
  dpi = 300
)

# ---------- GENERATE REPORT ----------
report_rmd <- file.path(reports_dir,
                        "Municipality_Spearman_Correlation_Heatmap_Report.Rmd")
)

cat("✓ Spearman correlation heatmap generated successfully\n")
cat("✓ Figure saved to:", figures_dir, "\n")
cat("✓ Report saved to:", reports_dir, "\n")
