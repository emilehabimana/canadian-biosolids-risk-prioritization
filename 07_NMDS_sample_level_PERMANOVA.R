# ============================================================
# PERMANOVA AND NMDS USING SPEARMAN CORRELATION DISTANCE
# ============================================================

# ---------- PATHS ----------
data_dir    <- "C:/Users/Admin/Documents/Data/Data"
tables_dir  <- "C:/Users/Admin/Documents/Data/Tables"
figures_dir <- "C:/Users/Admin/Documents/Data/Figures"
reports_dir <- "C:/Users/Admin/Documents/Data/Reports"

dir.create(tables_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(reports_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- LIBRARIES ----------
library(vegan)
library(ggplot2)
library(rmarkdown)

# ---------- LOAD DISTANCE MATRIX ----------
dist_df <- read.csv(
  file.path(data_dir,
            "SampleLevel_Spearman_Correlation_Distance_Matrix.csv"),
  row.names = 1,
  check.names = FALSE
)

# Clean empty column names
dist_df <- dist_df[, names(dist_df) != "" & !is.na(names(dist_df))]

dist_mat <- as.matrix(dist_df)
mode(dist_mat) <- "numeric"

# Enforce valid distance properties
dist_mat[!is.finite(dist_mat)] <- 0
dist_mat[dist_mat < 0] <- 0
dist_mat[dist_mat > 1] <- 1
dist_mat <- (dist_mat + t(dist_mat)) / 2
diag(dist_mat) <- 0

cor_dist <- as.dist(dist_mat)

# ---------- LOAD MUNICIPALITY METADATA ----------
ec <- read.csv(
  file.path(data_dir, "ec_concentration_logtransformed_PERMANOVA.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

municipality <- factor(toupper(trimws(ec[[1]])))

# Check alignment
if (length(municipality) != attr(cor_dist, "Size")) {
  stop("Metadata length does not match distance matrix size.")
}

# ---------- PERMANOVA ----------
set.seed(123)

perm <- adonis2(
  cor_dist ~ municipality,
  permutations = 999
)

perm_df <- as.data.frame(perm)

write.csv(
  perm_df,
  file.path(tables_dir,
            "PERMANOVA_Results_SpearmanDistance.csv")
)

# ---------- PERMDISP ----------
disp <- betadisper(cor_dist, municipality)
disp_test <- permutest(disp, permutations = 999)

disp_df <- data.frame(
  F_value = disp_test$tab[1, "F"],
  P_value = disp_test$tab[1, "Pr(>F)"]
)

write.csv(
  disp_df,
  file.path(tables_dir,
            "PERMDISP_Results_SpearmanDistance.csv"),
  row.names = FALSE
)

# ---------- NMDS ----------
set.seed(123)

nmds <- metaMDS(
  cor_dist,
  k = 2,
  trymax = 200,
  autotransform = FALSE
)

coords <- as.data.frame(scores(nmds))
coords$Municipality <- municipality

write.csv(
  coords,
  file.path(tables_dir,
            "NMDS_SampleLevel_Coordinates.csv"),
  row.names = FALSE
)

# ---------- NMDS FIGURE ----------
p <- ggplot(coords, aes(NMDS1, NMDS2, color = Municipality)) +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linewidth = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    color = "Municipality"
  )

ggsave(
  filename = "NMDS_SampleLevel_PERMANOVA_Spearman.png",
  plot = p,
  path = figures_dir,
  width = 10,
  height = 8,
  dpi = 300
)

# ---------- REPORT ----------
report_rmd <- file.path(reports_dir,
                        "PERMANOVA_SpearmanDistance_Report.Rmd")

report_text <- c(
  "---",
  "title: \"PERMANOVA of EC Co-occurrence Patterns by Municipality\"",
  "output: word_document",
  "---",
  "",
  "## Distance model",
  "",
  "Multivariate dissimilarity among samples was quantified using Spearman",
  "correlation distance, defined as d_ij = 1 - rho_ij, where rho_ij is the",
  "Spearman rank correlation between EC profiles of samples i and j.",
  "",
  "## PERMANOVA",
  "",
  "PERMANOVA was applied to test whether EC co-occurrence patterns differed",
  "among municipalities using 999 permutations.",
  "",
  "## PERMDISP",
  "",
  "Homogeneity of multivariate dispersion was assessed using PERMDISP.",
  "",
  "## Results",
  "",
  paste(
    "PERMANOVA results indicated",
    ifelse(perm_df$`Pr(>F)`[1] < 0.05,
           "statistically significant differences",
           "no statistically significant differences"),
    "in EC co-occurrence patterns among municipalities",
    "(pseudo-F =", round(perm_df$F[1], 2),
    ", R2 =", round(perm_df$R2[1], 3),
    ", p =", signif(perm_df$`Pr(>F)`[1], 3), ")."
  ),
  "",
  paste(
    "PERMDISP results showed",
    ifelse(disp_df$P_value > 0.05,
           "no significant differences in dispersion.",
           "significant differences in dispersion, requiring cautious interpretation.")
  ),
  "",
  "NMDS ordination based on the same distance matrix provided",
  "a visual representation of municipality-level separation.",
  "",
  "## References",
  "",
  "- Anderson, M.J. (2001). A new method for non-parametric multivariate analysis of variance. Austral Ecology.",
  "- Clarke, K.R. (1993). Non-parametric multivariate analyses of changes in community structure. Australian Journal of Ecology."
)

writeLines(report_text, report_rmd)

render(
  report_rmd,
  file.path(reports_dir,
            "PERMANOVA_SpearmanDistance_Report.docx"),
  quiet = TRUE
)

cat("✓ PERMANOVA, PERMDISP, NMDS, and report generated successfully\n")
