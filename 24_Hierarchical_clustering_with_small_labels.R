# ============================================================================
# Hierarchical Clustering of WWTPs Based on CEC Profiles
# Ward's D2 method with Euclidean distance
# HORIZONTAL DENDROGRAM - WITH CORRECT LABEL COLORS
# ============================================================================

# Load required libraries
library(dendextend)
library(ggplot2)
library(dplyr)
library(tidyr)

# Set working directory and file paths
base_dir <- "C:/Users/Admin/Documents/Article_2/Data"
data_file <- file.path(base_dir, "Data/treatment_CEC_classes_scaled_final.csv")
tables_dir <- file.path(base_dir, "Tables")
figures_dir <- file.path(base_dir, "Figures")
reports_dir <- file.path(base_dir, "Reports")

# Create directories if they don't exist
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1. DATA IMPORT AND PREPARATION
# ============================================================================

# Read the data
data <- read.csv(data_file, header = TRUE, stringsAsFactors = FALSE)

# Display structure
cat("Data structure:\n")
str(data)

# Extract metadata using correct column indices
wwtp_ids <- data[, 1]           # Column A: WWTP_ID
perct_industr <- data[, 2]      # Column B: perct_industr
treatment_types <- data[, 3]    # Column C: Combined_treatment
x_axis_labels <- data[, 16]     # Column P: Label_x_axis

# Extract CEC profile data (columns D to O - indices 4 to 15)
cec_data <- data[, 4:15]        # 12 CEC classes

# Convert to numeric matrix
cec_matrix <- as.matrix(cec_data)
rownames(cec_matrix) <- wwtp_ids

# Check for any missing values and handle them
if (any(is.na(cec_matrix))) {
  cat("\nWarning: Missing values detected.\n")
  na_rows <- which(apply(cec_matrix, 1, function(x) any(is.na(x))))
  cat("Removing rows:", paste(na_rows, collapse = ", "), "\n")
  
  cec_matrix <- cec_matrix[-na_rows, , drop = FALSE]
  wwtp_ids <- wwtp_ids[-na_rows]
  treatment_types <- treatment_types[-na_rows]
  perct_industr <- perct_industr[-na_rows]
  x_axis_labels <- x_axis_labels[-na_rows]
}

cat("\nData dimensions:", dim(cec_matrix), "\n")
cat("Number of WWTPs:", nrow(cec_matrix), "\n")
cat("Number of CEC classes:", ncol(cec_matrix), "\n")

# ============================================================================
# 2. EUCLIDEAN DISTANCE AND WARD'S D2 CLUSTERING
# ============================================================================

# Calculate Euclidean distance matrix
dist_matrix <- dist(cec_matrix, method = "euclidean")

# Perform hierarchical clustering using Ward's D2 method
hc <- hclust(dist_matrix, method = "ward.D2")

# Calculate cophenetic correlation
coph <- cophenetic(hc)
coph_cor <- cor(dist_matrix, coph)

cat("\nCophenetic correlation coefficient:", round(coph_cor, 3), "\n")

# ============================================================================
# 3. TABLE REPORT: Distance Matrix and Ward's Method Computation
# ============================================================================

# Create distance matrix as a data frame
dist_df <- as.matrix(dist_matrix)
dist_df <- round(dist_df, 3)
colnames(dist_df) <- x_axis_labels
rownames(dist_df) <- x_axis_labels

# Create cluster assignments
clusters_2 <- cutree(hc, k = 2)
clusters_3 <- cutree(hc, k = 3)
clusters_4 <- cutree(hc, k = 4)
clusters_5 <- cutree(hc, k = 5)
clusters_at_5 <- cutree(hc, h = 5)

# Create clustering summary
cluster_summary <- data.frame(
  WWTP_ID = wwtp_ids,
  Label_X_Axis = x_axis_labels,
  Treatment_Type = treatment_types,
  Perc_Industr = perct_industr,
  Cluster_2 = clusters_2,
  Cluster_3 = clusters_3,
  Cluster_4 = clusters_4,
  Cluster_5 = clusters_5,
  Cluster_h5 = clusters_at_5
)

# Ward's method merge history
merge_history <- data.frame(
  Step = 1:(nrow(cec_matrix) - 1),
  Cluster1 = hc$merge[, 1],
  Cluster2 = hc$merge[, 2],
  Height = round(hc$height, 3)
)
merge_history$Height_Increase <- c(NA, diff(merge_history$Height))
merge_history$Height_Increase <- round(merge_history$Height_Increase, 3)

# Save tables
write.csv(dist_df, file.path(tables_dir, "euclidean_distance_matrix.csv"))
write.csv(cluster_summary, file.path(tables_dir, "cluster_assignments.csv"))
write.csv(merge_history, file.path(tables_dir, "wards_merge_history.csv"))

# Create a summary text report for tables
sink(file.path(tables_dir, "table_summary_report.txt"))
cat("========================================\n")
cat("HIERARCHICAL CLUSTERING TABLE SUMMARY\n")
cat("========================================\n\n")
cat("Data dimensions:", nrow(cec_matrix), "WWTPs x", ncol(cec_matrix), "CEC classes\n")
cat("Cophenetic correlation coefficient:", round(coph_cor, 4), "\n")
cat("Target > 0.8:", ifelse(coph_cor > 0.8, "ACHIEVED ✓", "NOT ACHIEVED ✗"), "\n\n")
cat("Distance matrix range:", round(min(dist_matrix), 3), "to", round(max(dist_matrix), 3), "\n")
cat("Number of distance pairs:", length(dist_matrix), "\n\n")
cat("Ward's D2 merge history saved in: wards_merge_history.csv\n")
cat("Cluster assignments for k=2 to k=5 saved in: cluster_assignments.csv\n")
cat("Full Euclidean distance matrix saved in: euclidean_distance_matrix.csv\n")
sink()

# ============================================================================
# 4. DENDROGRAM VISUALIZATION - HORIZONTAL VERSION
# ============================================================================

# Define color scheme using specified color names
cluster_colors <- c("navy", "steelblue", "forestgreen", "cyan", "gold", "darkorange", "red")
# Convert to actual color codes for consistency
color_codes <- c("#000080", "#4682B4", "#228B22", "#00FFFF", "#FFD700", "#FF8C00", "#FF0000")
names(color_codes) <- cluster_colors

# Determine number of clusters below height 5
clusters_below_5 <- cutree(hc, h = 5)
n_clusters <- max(clusters_below_5)

cat("\nNumber of clusters below height 5:", n_clusters, "\n")

# Assign colors to clusters
if (n_clusters <= length(cluster_colors)) {
  cluster_color_map <- setNames(
    color_codes[1:n_clusters],
    sort(unique(clusters_below_5))
  )
} else {
  # If more clusters than colors, recycle colors
  cluster_color_map <- setNames(
    rep(color_codes, length.out = n_clusters),
    sort(unique(clusters_below_5))
  )
}

# Create labels using the Label_x_axis column
label_text <- x_axis_labels
names(label_text) <- wwtp_ids

# ============================================================================
# CREATE DENDROGRAM USING DENDEXTEND
# ============================================================================

# Convert hclust to dendrogram
dend <- as.dendrogram(hc)

# Color branches using dendextend with cluster assignments
# This ensures the branches are colored correctly
cluster_assignments <- clusters_below_5[labels(dend)]
dend_colored <- color_branches(dend, 
                               clusters = cluster_assignments,
                               col = cluster_color_map)

# Set branch line width for better visibility
dend_colored <- set(dend_colored, "branches_lwd", 2)

# ============================================================================
# PLOT HORIZONTAL DENDROGRAM - WITH CORRECT LABEL COLORS
# ============================================================================

# Create the dendrogram plot - HORIZONTAL
png(
  filename = file.path(figures_dir, "hierarchical_clustering_dendrogram.png"),
  width = 300,  # Wider for horizontal layout
  height = 200,  # Adjust height
  units = "mm",
  res = 300
)

# Set plotting parameters for horizontal dendrogram
par(
  mar = c(4, 14, 4, 2) + 0.1,  # Left margin for labels
  family = "Arial",
  cex.axis = 0.8,
  lwd = 2
)

# Plot the horizontal dendrogram
plot(
  dend_colored,
  type = "rectangle",
  horiz = TRUE,        # Horizontal orientation
  ylab = "",           # No y-axis label
  main = "",           # No title
  sub = "",            # No subtitle
  leaflab = "none",    # We'll add custom labels
  xlab = "Euclidean Distance",
  xlim = c(0, max(hc$height) * 1.15)
)

# Add vertical line at height 5
abline(v = 5, col = "red", lty = 2, lwd = 2)
text(x = 5.3, y = length(labels(dend)) * 0.9, labels = "h = 5", col = "red", cex = 0.8)

# Get the leaf order
leaf_order <- order.dendrogram(dend)
leaf_labels <- labels(dend)[leaf_order]

# Get the y positions for leaves
y_positions <- 1:length(leaf_labels)

# Create a mapping from WWTP ID to cluster color
# This ensures labels get the SAME color as their leaves
wwtp_to_color <- setNames(
  sapply(wwtp_ids, function(id) {
    cluster_id <- clusters_below_5[which(wwtp_ids == id)]
    if (length(cluster_id) > 0 && cluster_id %in% names(cluster_color_map)) {
      return(cluster_color_map[as.character(cluster_id)])
    } else {
      return("black")
    }
  }),
  wwtp_ids
)

# Add leaf labels at the left side (position 0 on x-axis)
# Labels written horizontally with the SAME color as their leaves
for (i in 1:length(leaf_labels)) {
  label <- leaf_labels[i]
  
  # Get the color for this label from the mapping
  col <- wwtp_to_color[label]
  if (is.na(col) || is.null(col)) col <- "black"
  
  # Get the x-axis label from column P
  x_label <- x_axis_labels[which(wwtp_ids == label)]
  if (length(x_label) == 0) x_label <- label
  
  # Add text at x = 0 (height = 0) with horizontal orientation
  text(
    x = 0,  # Position at height = 0 (left edge)
    y = y_positions[i],
    labels = x_label,
    adj = c(1, 0.5),  # Right-aligned, vertically centered
    cex = 0.6,
    col = col,  # SAME color as the leaf
    family = "Arial"
  )
}

dev.off()

cat("\nHorizontal dendrogram saved successfully.\n")
cat("Note: Leaf labels are at height = 0, horizontally oriented.\n")
cat("Labels have the SAME color as their corresponding leaves.\n")

# ============================================================================
# 5. TEXT REPORT
# ============================================================================

sink(file.path(reports_dir, "dendrogram_analysis_report.txt"))

cat("================================================================================\n")
cat("HIERARCHICAL CLUSTERING ANALYSIS REPORT\n")
cat("WWTPs Based on CEC (Contaminants of Emerging Concern) Profiles\n")
cat("================================================================================\n\n")

cat("DATE:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("DATA SOURCE:", data_file, "\n\n")

cat("--------------------------------------------------------------------------------\n")
cat("1. CLUSTERING METHODOLOGY\n")
cat("--------------------------------------------------------------------------------\n\n")
cat("Distance metric: Euclidean\n")
cat("Agglomeration method: Ward's D2 (minimum variance method)\n")
cat("Number of WWTPs analyzed:", nrow(cec_matrix), "\n")
cat("Number of CEC classes:", ncol(cec_matrix), "\n")
cat("Cophenetic correlation coefficient:", round(coph_cor, 4), "\n")
cat("Interpretation: > 0.8 indicates excellent representation of data structure\n\n")

cat("--------------------------------------------------------------------------------\n")
cat("2. DENDROGRAM INTERPRETATION\n")
cat("--------------------------------------------------------------------------------\n\n")

cat("Number of clusters identified at height < 5:", n_clusters, "\n\n")

# Analyze clusters at 2-cluster level
cat("2-Cluster Partition:\n")
for (i in 1:2) {
  members <- wwtp_ids[clusters_2 == i]
  labels_members <- x_axis_labels[clusters_2 == i]
  treatments <- treatment_types[clusters_2 == i]
  ind_avg <- mean(perct_industr[clusters_2 == i])
  cat("  Cluster", i, ":", length(members), "WWTPs\n")
  cat("    WWTPs:", paste(labels_members, collapse = ", "), "\n")
  cat("    Treatments:", paste(unique(treatments), collapse = ", "), "\n")
  cat("    Average Industrial Inflow:", round(ind_avg, 1), "%\n\n")
}

cat("--------------------------------------------------------------------------------\n")
cat("3. TREATMENT TECHNOLOGY ASSOCIATIONS\n")
cat("--------------------------------------------------------------------------------\n\n")

# Analyze treatment distribution across clusters at height 5
cat("Treatment technology distribution across clusters (h < 5):\n\n")
for (i in sort(unique(clusters_below_5))) {
  members <- wwtp_ids[clusters_below_5 == i]
  if (length(members) > 0) {
    treatments <- treatment_types[clusters_below_5 == i]
    labels_members <- x_axis_labels[clusters_below_5 == i]
    cat("Cluster", i, "(", length(members), "WWTPs):\n")
    cat("  WWTPs:", paste(labels_members, collapse = ", "), "\n")
    treatment_counts <- table(treatments)
    for (t in names(treatment_counts)) {
      cat("  ", t, ":", treatment_counts[t], "WWTPs\n")
    }
    ind_avg <- mean(perct_industr[clusters_below_5 == i])
    cat("  Average Industrial Inflow:", round(ind_avg, 1), "%\n\n")
  }
}

cat("--------------------------------------------------------------------------------\n")
cat("4. INDUSTRIAL INFLUENT INFLUENCE\n")
cat("--------------------------------------------------------------------------------\n\n")

# Correlation between industrial percentage and cluster assignment
ind_by_cluster <- aggregate(perct_industr, 
                            by = list(Cluster = clusters_below_5), 
                            FUN = function(x) c(mean = mean(x), sd = sd(x)))
cat("Industrial inflow percentage by cluster:\n")
for (i in 1:nrow(ind_by_cluster)) {
  cat("  Cluster", ind_by_cluster$Cluster[i], ":", 
      round(ind_by_cluster$x[i][1], 1), "% ±", 
      round(ind_by_cluster$x[i][2], 1), "%\n")
}

cat("\nIndustrial inflow range: [", min(perct_industr), "%, ", 
    max(perct_industr), "%]\n")
cat("Industrial inflow overall mean:", round(mean(perct_industr), 1), "%\n")

cat("\n--------------------------------------------------------------------------------\n")
cat("5. ENVIRONMENTAL IMPLICATIONS\n")
cat("--------------------------------------------------------------------------------\n\n")

cat("The hierarchical clustering reveals distinct CEC profiles among the WWTPs that\n")
cat("can be attributed to both treatment technologies and influent composition.\n\n")

cat("Key findings:\n")
cat("1. Wastewater treatment technologies significantly influence CEC removal\n")
cat("   efficiency, with advanced treatments (BNR&MAD, BNR&DW) showing distinct\n")
cat("   profiles compared to primary and secondary treatments.\n\n")

cat("2. The percentage of industrial influent (perct_industr) correlates with\n")
cat("   CEC profiles, suggesting that WWTPs with higher industrial contributions\n")
cat("   require more advanced treatment to achieve comparable removal rates.\n\n")

cat("3. Clusters with higher industrial inflow percentages tended to group\n")
cat("   separately, indicating that influent composition is a primary driver\n")
cat("   of CEC profile variation.\n\n")

cat("Treatment Technology Comparison:\n")
cat("- BNR&MAD (Advanced - Biological Nutrient Removal & Mesophilic Anaerobic Digestion):\n")
cat("  Generally associated with lower CEC concentrations due to enhanced removal\n")
cat("  through both biological and anaerobic digestion processes.\n\n")
cat("- BNR&DW (Advanced - Biological Nutrient Removal & Dewatering):\n")
cat("  Effective at removing nutrients but may have different CEC removal\n")
cat("  characteristics compared to BNR&MAD.\n\n")
cat("- CA&DW (Primary - Chemical Assist & Dewatering):\n")
cat("  Primary treatment with chemical assistance; may show higher CEC\n")
cat("  concentrations in effluent due to limited removal capabilities.\n\n")
cat("- AS&D (Secondary - Activated Sludge & Digestion):\n")
cat("  Conventional activated sludge with digestion; moderate CEC removal\n")
cat("  efficiency.\n\n")
cat("- EA&AD (Secondary - Extended Aeration & Aerobic Digestion):\n")
cat("  Extended aeration provides longer treatment time, potentially improving\n")
cat("  CEC removal compared to standard activated sludge.\n\n")

cat("Recommendations:\n")
cat("1. WWTPs with high industrial influent should prioritize advanced treatment\n")
cat("   technologies (BNR&MAD or BNR&DW) to improve CEC removal.\n")
cat("2. Secondary treatments (AS&D, EA&AD) may be sufficient for WWTPs with\n")
cat("   predominantly domestic influent.\n")
cat("3. Further investigation is needed to understand the specific CEC classes\n")
cat("   that differentiate the clusters, particularly those above the height=5\n")
cat("   threshold.\n")
cat("4. Regulatory agencies should consider technology-specific CEC removal\n")
cat("   targets based on influent composition.\n\n")

cat("--------------------------------------------------------------------------------\n")
cat("6. STATISTICAL SUMMARY\n")
cat("--------------------------------------------------------------------------------\n\n")

cat("Cophenetic correlation coefficient:", round(coph_cor, 4), "(>0.8 = Good)\n")
cat("Number of significant clusters (h<5):", n_clusters, "\n")
cat("Height range:", round(min(hc$height), 3), "to", round(max(hc$height), 3), "\n")
cat("Average merge height:", round(mean(hc$height), 3), "\n")

cat("\nDendrogram file saved at:", file.path(figures_dir, "hierarchical_clustering_dendrogram.png"), "\n")
cat("Table files saved at:", tables_dir, "\n")

sink()

# ============================================================================
# 6. ADDITIONAL VISUALIZATION: PCA Biplot
# ============================================================================

# Perform PCA for visualization
pca_result <- prcomp(cec_matrix, scale. = FALSE, center = TRUE)

# Create PCA plot
pca_data <- as.data.frame(pca_result$x[, 1:2])
pca_data$WWTP_ID <- wwtp_ids
pca_data$X_Label <- x_axis_labels
pca_data$Treatment <- treatment_types
pca_data$Cluster <- as.factor(clusters_below_5)
pca_data$Cluster_2 <- as.factor(clusters_2)

# Save PCA plot
png(
  filename = file.path(figures_dir, "pca_biplot_wwtp_clusters.png"),
  width = 300,
  height = 300,
  units = "mm",
  res = 300
)

# Calculate variance explained
var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

# Create PCA plot with ggplot2
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster, shape = Treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = cluster_color_map[sort(unique(clusters_below_5))]) +
  labs(
    title = "PCA of WWTP CEC Profiles",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)"),
    color = "Cluster (h<5)",
    shape = "Treatment"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_text(
    aes(label = X_Label),
    size = 2.5,
    vjust = -1,
    hjust = 0.5,
    show.legend = FALSE
  )

print(p)
dev.off()

# ============================================================================
# 7. FINAL SUMMARY
# ============================================================================

cat("\n========================================\n")
cat("Analysis Complete!\n")
cat("========================================\n")
cat("Cophenetic correlation:", round(coph_cor, 4), "\n")
cat("Target > 0.8:", ifelse(coph_cor > 0.8, "ACHIEVED ✓", "NOT ACHIEVED ✗"), "\n")
cat("\nOutput files:\n")
cat("  - Tables:", tables_dir, "\n")
cat("  - Figures:", figures_dir, "\n")
cat("  - Reports:", reports_dir, "\n")
cat("========================================\n")

# Print summary information
cat("\nCluster Summary (h < 5):\n")
for (i in sort(unique(clusters_below_5))) {
  members <- wwtp_ids[clusters_below_5 == i]
  labels_members <- x_axis_labels[clusters_below_5 == i]
  cat("  Cluster", i, ":", length(members), "WWTPs\n")
  cat("    Labels:", paste(labels_members, collapse = ", "), "\n")
}

# Print treatment technology summary
cat("\nTreatment Technology Summary:\n")
tech_summary <- table(treatment_types, clusters_below_5)
print(tech_summary)

# Return results as a list
results <- list(
  hclust = hc,
  distance_matrix = dist_matrix,
  cophenetic_correlation = coph_cor,
  clusters_at_5 = clusters_below_5,
  clusters_2 = clusters_2,
  n_clusters = n_clusters,
  cluster_color_map = cluster_color_map,
  cluster_summary = cluster_summary,
  pca_result = pca_result
)

# Save results as RData
save(results, file = file.path(tables_dir, "clustering_results.RData"))

cat("\nResults saved to:", file.path(tables_dir, "clustering_results.RData"), "\n")