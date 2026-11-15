# ============================================================
# 3D NMDS WITH MUNICIPALITY CLUSTERS AND EQUAL AXIS SCALING
# ============================================================

# ---------- PATHS ----------
data_dir    <- "C:/Users/Admin/Documents/Article_2/Data/Data"
figures_dir <- "C:/Users/Admin/Documents/Article_2/Data/Figures"

dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- LOAD DISTANCE MATRIX ----------
dist_df <- read.csv(
  file.path(
    data_dir,
    "SampleLevel_Spearman_Correlation_Distance_Matrix.csv"
  ),
  row.names = 1,
  check.names = FALSE
)

dist_mat <- as.matrix(dist_df)
mode(dist_mat) <- "numeric"

# Ensure valid distance properties
dist_mat[!is.finite(dist_mat)] <- 0
dist_mat[dist_mat < 0] <- 0
dist_mat[dist_mat > 1] <- 1
dist_mat <- (dist_mat + t(dist_mat)) / 2
diag(dist_mat) <- 0

cor_dist <- as.dist(dist_mat)

# ---------- LOAD MUNICIPALITY METADATA ----------
meta <- read.csv(
  file.path(
    data_dir,
    "ec_concentration_logtransformed_PERMANOVA.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

municipality <- factor(toupper(trimws(meta[[1]])))

# ---------- 3D NMDS ----------
set.seed(123)

nmds_3d <- metaMDS(
  cor_dist,
  k = 3,
  trymax = 300,
  autotransform = FALSE
)

coords <- as.data.frame(scores(nmds_3d))
coords$Municipality <- municipality

# ---------- AUTOMATIC MUNICIPALITY CLUSTERING ----------
hc <- hclust(cor_dist, method = "average")
coords$Cluster <- factor(cutree(hc, k = 4))

# ---------- EQUAL AXIS SCALING ----------
ranges <- apply(coords[, 1:3], 2, range)
max_range <- max(ranges[2, ] - ranges[1, ])
centers <- colMeans(ranges)

axis_limits <- rbind(
  centers - max_range / 2,
  centers + max_range / 2
)

# ---------- INTERACTIVE 3D NMDS WITH CLUSTERS ----------
p3d <- plot_ly(
  coords,
  x = ~NMDS1,
  y = ~NMDS2,
  z = ~NMDS3,
  color = ~Cluster,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 4)
) %>%
  layout(
    scene = list(
      xaxis = list(
        title = "NMDS1",
        range = axis_limits[, 1]
      ),
      yaxis = list(
        title = "NMDS2",
        range = axis_limits[, 2]
      ),
      zaxis = list(
        title = "NMDS3",
        range = axis_limits[, 3]
      )
    )
  )

htmlwidgets::saveWidget(
  p3d,
  file.path(
    figures_dir,
    "NMDS_3D_Municipality_Clusters_interactive.html"
  )
)

# ---------- STATIC 3D NMDS WITH CLUSTERS ----------
cluster_colors <- rainbow(length(levels(coords$Cluster)))

png(
  file.path(
    figures_dir,
    "NMDS_3D_Municipality_Clusters_static.png"
  ),
  width = 1200,
  height = 1000,
  res = 150
)

scatterplot3d(
  coords$NMDS1,
  coords$NMDS2,
  coords$NMDS3,
  color = cluster_colors[coords$Cluster],
  pch = 19,
  xlab = "NMDS1",
  ylab = "NMDS2",
  zlab = "NMDS3",
  xlim = axis_limits[, 1],
  ylim = axis_limits[, 2],
  zlim = axis_limits[, 3]
)

legend(
  "topright",
  legend = paste("Cluster", levels(coords$Cluster)),
  col = cluster_colors,
  pch = 19,
  cex = 0.9
)

dev.off()

cat("✓ 3D NMDS with clusters and equal axis scaling completed\n")
cat("✓ Stress:", round(nmds_3d$stress, 3), "\n")
``