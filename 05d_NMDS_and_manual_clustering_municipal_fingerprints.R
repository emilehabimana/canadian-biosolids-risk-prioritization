# ============================================================
# NMDS OF MUNICIPAL EC FINGERPRINTS
# AUTOMATIC CLUSTERING INTO 4 CLUSTERS + ELLIPSES
# ============================================================

# ---------- PATHS ----------
data_dir    <- "C:/Users/Admin/Documents/Article_2/Data/Data"
tables_dir  <- "C:/Users/Admin/Documents/Article_2/Data/Tables"
figures_dir <- "C:/Users/Admin/Documents/Article_2/Data/Figures"

dir.create(tables_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- LIBRARIES ----------
library(vegan)
library(ggplot2)

# ---------- LOAD SPEARMAN CORRELATION DISTANCE MATRIX ----------
dist_df <- read.csv(
  file.path(
    data_dir,
    "Municipality_Spearman_Correlation_Distance_MUN.csv"
  ),
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dist_mat <- as.matrix(dist_df)
mode(dist_mat) <- "numeric"

# ---------- SANITIZE DISTANCE MATRIX ----------
dist_mat[!is.finite(dist_mat)] <- 0
dist_mat[dist_mat < 0] <- 0
dist_mat[dist_mat > 1] <- 1
dist_mat <- (dist_mat + t(dist_mat)) / 2
diag(dist_mat) <- 0

cor_dist <- as.dist(dist_mat)

# ---------- NMDS ----------
set.seed(123)

nmds <- metaMDS(
  cor_dist,
  k = 2,
  trymax = 200,
  autotransform = FALSE
)

coords <- as.data.frame(scores(nmds))
coords$Municipality <- rownames(coords)

# ---------- AUTOMATIC HIERARCHICAL CLUSTERING ----------
hc <- hclust(cor_dist, method = "average")

# Number of clusters
k <- 4
coords$Cluster <- factor(cutree(hc, k = k))

# ---------- SAVE CLUSTER ASSIGNMENTS ----------
write.csv(
  coords[, c("Municipality", "Cluster")],
  file.path(
    tables_dir,
    "Municipality_EC_AutoCluster_Assignments.csv"
  ),
  row.names = FALSE
)

# ---------- NMDS PLOT WITH ELLIPSES ----------
p <- ggplot(
  coords,
  aes(NMDS1, NMDS2, color = Cluster)
) +
  geom_point(size = 4) +
  geom_text(
    aes(label = Municipality),
    size = 5,
    fontface = "bold",
    vjust = -0.8,
    show.legend = FALSE
  ) +
  stat_ellipse(
    aes(group = Cluster),
    type = "t",
    linetype = "dashed",
    linewidth = 1.2,
    show.legend = FALSE
  ) +
  theme_bw() +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    color = "Cluster"
  ) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    
    # REMOVE GRID LINES
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.title = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.text = element_text(
      size = 15,
      face = "bold",
      color = "black"
    ),
    
    legend.title = element_text(
      size = 16,
      face = "bold"
    ),
    legend.text = element_text(
      size = 16
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 6)
    )
  )

# ---------- SAVE FIGURE ----------
ggsave(
  filename = "NMDS_Municipality_EC_Fingerprints_AutoClusters_Ellipses.png",
  plot = p,
  path = figures_dir,
  width = 12,
  height = 10,
  dpi = 300
)

cat("✓ NMDS with automatic clustering (k = 4) and ellipses completed successfully\n")
``