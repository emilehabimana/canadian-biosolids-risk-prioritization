# ============================================================
# 3D NMDS WITH CLUSTER ELLIPSOIDS AND TRAJECTORIES (SAFE)
# ============================================================

# ---------- PATHS ----------
data_dir <- "C:/Users/Admin/Documents/Data/Data"
figures_dir <- "C:/Users/Admin/Documents/Data/Figures"

if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

# ---------- LIBRARIES ----------
library(vegan)
library(plotly)

# ---------- LOAD DISTANCE MATRIX ----------
dist_df <- read.csv(
  file.path(data_dir,
            "SampleLevel_Spearman_Correlation_Distance_Matrix.csv"),
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dist_df <- dist_df[, names(dist_df) != "" & !is.na(names(dist_df))]
dist_mat <- as.matrix(dist_df)
mode(dist_mat) <- "numeric"

dist_mat[!is.finite(dist_mat)] <- 0
dist_mat[dist_mat < 0] <- 0
dist_mat[dist_mat > 1] <- 1
dist_mat <- (dist_mat + t(dist_mat)) / 2
diag(dist_mat) <- 0

cor_dist <- as.dist(dist_mat)

# ---------- LOAD METADATA ----------
meta <- read.csv(
  file.path(data_dir,
            "ec_concentration_logtransformed_PERMANOVA.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

municipality <- factor(toupper(trimws(meta[[1]])))

if (length(municipality) != attr(cor_dist, "Size")) {
  stop("Mismatch between distance matrix and metadata length.")
}

# ---------- 3D NMDS ----------
set.seed(123)
nmds3 <- metaMDS(cor_dist, k = 3, trymax = 300, autotransform = FALSE)

coords <- as.data.frame(scores(nmds3))
coords$Cluster <- factor(cutree(hclust(cor_dist, method = "average"), k = 4))

# ---------- EQUAL AXIS SCALING ----------
xyz <- as.matrix(coords[, c("NMDS1", "NMDS2", "NMDS3")])
ranges <- apply(xyz, 2, range)
span <- max(ranges[2, ] - ranges[1, ])
center <- colMeans(ranges)

xlim <- center[1] + c(-span, span) / 2
ylim <- center[2] + c(-span, span) / 2
zlim <- center[3] + c(-span, span) / 2

# ---------- FUNCTION TO BUILD ELLIPSOID ----------
make_ellipsoid <- function(mu, Sigma, level = 0.95, n = 40) {
  
  eig <- eigen(Sigma)
  radii <- sqrt(qchisq(level, df = 3)) * sqrt(eig$values)
  
  u <- seq(0, 2 * pi, length.out = n)
  v <- seq(0, pi, length.out = n)
  
  xyz_array <- array(0, dim = c(n, n, 3))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      xyz_array[i, j, ] <- c(
        radii[1] * cos(u[i]) * sin(v[j]),
        radii[2] * sin(u[i]) * sin(v[j]),
        radii[3] * cos(v[j])
      )
    }
  }
  
  rot <- eig$vectors
  for (i in 1:n) {
    for (j in 1:n) {
      xyz_array[i, j, ] <- rot %*% xyz_array[i, j, ] + mu
    }
  }
  
  list(
    x = xyz_array[, , 1],
    y = xyz_array[, , 2],
    z = xyz_array[, , 3]
  )
}

# ---------- BUILD 3D NMDS PLOT ----------
p <- plot_ly()

colors <- c("red", "green", "blue", "purple")

for (k in levels(coords$Cluster)) {
  
  idx <- coords$Cluster == k
  pts <- xyz[idx, , drop = FALSE]
  
  # Scatter points
  p <- add_trace(
    p,
    x = pts[, 1],
    y = pts[, 2],
    z = pts[, 3],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 4),
    name = paste("Cluster", k),
    color = I(colors[as.numeric(k)])
  )
  
  # Centroid
  mu <- colMeans(pts)
  p <- add_trace(
    p,
    x = mu[1],
    y = mu[2],
    z = mu[3],
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 6, symbol = "x"),
    showlegend = FALSE,
    color = I(colors[as.numeric(k)])
  )
  
  # Trajectories
  for (i in seq_len(nrow(pts))) {
    p <- add_trace(
      p,
      x = c(mu[1], pts[i, 1]),
      y = c(mu[2], pts[i, 2]),
      z = c(mu[3], pts[i, 3]),
      type = "scatter3d",
      mode = "lines",
      line = list(width = 1),
      showlegend = FALSE,
      color = I(colors[as.numeric(k)])
    )
  }
  
  # Ellipsoid (if enough points)
  if (nrow(pts) >= 4) {
    ell <- make_ellipsoid(mu, cov(pts))
    p <- add_surface(
      p,
      x = ell$x,
      y = ell$y,
      z = ell$z,
      opacity = 0.15,
      showscale = FALSE,
      showlegend = FALSE,
      surfacecolor = matrix(as.numeric(k), nrow(ell$x)),
      colorscale = list(c(0, 1), c(colors[as.numeric(k)], colors[as.numeric(k)]))
    )
  }
}

# ---------- LAYOUT WITH ARIAL BLACK FONT FOR ALL AXIS LABELS ----------
p <- layout(
  p,
  scene = list(
    xaxis = list(
      title = list(
        text = "NMDS1",
        font = list(family = "Arial", size = 14, color = "black")
      ),
      tickfont = list(family = "Arial", size = 12, color = "black"),
      range = xlim,
      color = "black"
    ),
    yaxis = list(
      title = list(
        text = "NMDS2",
        font = list(family = "Arial", size = 14, color = "black")
      ),
      tickfont = list(family = "Arial", size = 12, color = "black"),
      range = ylim,
      color = "black"
    ),
    zaxis = list(
      title = list(
        text = "NMDS3",
        font = list(family = "Arial", size = 14, color = "black")
      ),
      tickfont = list(family = "Arial", size = 12, color = "black"),
      range = zlim,
      color = "black"
    )
  ),
  # Global font settings for any additional text elements
  font = list(family = "Arial", color = "black")
)

# ---------- SAVE INTERACTIVE FIGURE ----------
htmlwidgets::saveWidget(
  p,
  file.path(figures_dir, "NMDS_3D_clusters_ellipsoids_SAFE.html")
)

cat("✓ 3D NMDS with ellipsoids and trajectories completed\n")
cat("✓ Stress:", round(nmds3$stress, 3), "\n")
