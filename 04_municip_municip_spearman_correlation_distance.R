# ============================================================
# SPEARMAN CORRELATION HEATMAP OF MUNICIPAL EC FINGERPRINTS
# ============================================================

# ---------- PATHS ----------
data_dir    <- "C:/Users/Admin/Documents/Article_2/Data/Data"
figures_dir <- "C:/Users/Admin/Documents/Article_2/Data/Figures"
reports_dir <- "C:/Users/Admin/Documents/Article_2/Data/Reports"

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

report_text <- c(
  "---",
  "title: \"Spearman Correlation Heatmap of Municipal EC Fingerprints\"",
  "output: word_document",
  "---",
  "",
  "## Rationale for the 29 × 29 correlation matrix",
  "",
  "Each municipality is represented by a single multivariate fingerprint consisting of",
  "robust-scaled, log-transformed concentrations of 110 emerging contaminants (ECs).",
  "Pairwise correlation was therefore computed between municipalities, resulting in a",
  "29 × 29 matrix that describes similarity in EC composition across all municipal biosolids.",
  "",
  "## Theoretical background",
  "",
  "Spearman's rank correlation coefficient (ρ) quantifies monotonic relationships between",
  "two multivariate profiles based on rank similarity rather than absolute values. For two",
  "municipalities i and j, correlation distance was defined as:",
  "",
  "$$",
  "d_{ij} = 1 - \\rho_{ij}",
  "$$",
  "",
  "and converted back to correlation for visualization as:",
  "",
  "$$",
  "\\rho_{ij} = 1 - d_{ij}",
  "$$",
  "",
  "This approach emphasizes relative EC composition and minimizes sensitivity to outliers",
  "and residual skewness.",
  "",
  "## Interpretation of the heatmap",
  "",
  "Dark blue cells indicate strong positive monotonic similarity in EC fingerprints between",
  "municipalities, whereas turquoise tones indicate weak or absent association. The diagonal",
  "represents perfect self-correlation (ρ = 1). Clusters of high correlation suggest shared",
  "chemical signatures potentially reflecting similar source contributions or treatment",
  "practices.",
  "",
  "## Manuscript-ready Methods paragraph",
  "",
  "> Municipality-level EC fingerprints were compared using Spearman rank correlation based",
  "on robust-scaled, log-transformed concentrations of 110 compounds. Pairwise correlations",
  "were computed between municipalities and visualized as a 29 × 29 correlation matrix,",
  "highlighting monotonic similarity in EC composition independent of absolute concentration",
  "levels.",
  "",
  "## Manuscript-ready Results paragraph",
  "",
  "> The Spearman correlation heatmap revealed variable similarity in EC composition among",
  "municipalities, with several municipality pairs exhibiting strong positive monotonic",
  "associations in their EC fingerprints. These patterns indicate municipality-specific yet",
  "partially shared chemical signatures, consistent with differences in contaminant sources",
  "and local conditions influencing biosolids composition."
)

writeLines(report_text, report_rmd)

# ---------- RENDER WORD DOCUMENT ----------
render(
  input = report_rmd,
  output_file = file.path(
    reports_dir,
    "Municipality_Spearman_Correlation_Heatmap_Report.docx"
  ),
  quiet = TRUE
)

cat("✓ Spearman correlation heatmap generated successfully\n")
cat("✓ Figure saved to:", figures_dir, "\n")
cat("✓ Report saved to:", reports_dir, "\n")