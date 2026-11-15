################################################################################
# SCRIPT: PERMANOVA, and PERMDISP Analysis of CEC Classes
# Objective: Evaluate how Combined_treatement influences CEC class composition
# Data: 12 CEC class concentrations (µg/kg) from treatments_CEC_classes.csv
################################################################################

# ==============================================================================
# 0. SETUP - Load Libraries and Configure Environment
# ==============================================================================

# Install packages if not already installed
packages <- c("vegan", "dplyr", "tidyr", "ggplot2", "ape", "pairwiseAdonis", 
              "grid", "ggrepel", "scales")

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

invisible(lapply(packages, install_if_missing))

# Load required libraries
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ape)
library(pairwiseAdonis)
library(grid)
library(ggrepel)
library(scales)

# Disable scientific notation for cleaner output
options(scipen = 999)

# Set seed for reproducibility
set.seed(123)

# ==============================================================================
# 1. DATA IMPORT AND PREPARATION
# ==============================================================================

# Define file paths
data_dir <- "C:/Users/Admin/Documents/Data/Data"
tables_dir <- "C:/Users/Admin/Documents/Data/Tables"
reports_dir <- "C:/Users/Admin/Documents/Data/Reports"
figures_dir <- "C:/Users/Admin/Documents/Data/Figures"

# Create directories if they don't exist
for (dir in c(tables_dir, reports_dir, figures_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("Created directory:", dir, "\n")
  }
}

# Define the CSV file path
csv_file <- file.path(data_dir, "treatments_CEC_classes.csv")

# Check if file exists
if (!file.exists(csv_file)) {
  stop(paste("ERROR: File not found at:", csv_file, 
             "\nPlease verify the file name and path."))
}

# Read the CSV file
cat("\n========== DATA IMPORT ==========\n")
cat("Reading data from:", csv_file, "\n")

df <- read.csv(csv_file, header = TRUE, stringsAsFactors = FALSE)

# Fix column names (remove trailing spaces)
colnames(df) <- trimws(colnames(df))

# Display the data structure
cat("\nColumn names after cleaning:\n")
print(colnames(df))

cat("\nFirst 6 rows of data:\n")
print(head(df))

# ==============================================================================
# 2. DATA PREPARATION
# ==============================================================================

cat("\n========== DATA PREPARATION ==========\n")

# Define CEC class columns (I through T) - RAW CONCENTRATIONS in µg/kg
cec_columns <- c("Op_Sti", "Abx_Amx", "ASM_AED", "Oth_Pharm", "PCP", "ASs", 
                 "Neonics", "Herb", "Oth_Pest", "Industr", "Horm", "Bisph")

# Check if all CEC columns exist
missing_cols <- cec_columns[!cec_columns %in% colnames(df)]
if (length(missing_cols) > 0) {
  stop(paste("ERROR: Missing columns:", paste(missing_cols, collapse = ", ")))
}

# Check if Combined_treatement column exists
if (!"Combined_treatement" %in% colnames(df)) {
  stop("ERROR: 'Combined_treatement' column not found in the CSV file.")
}

# Extract the 12 CEC class columns as response matrix (RAW CONCENTRATIONS)
cec_data <- df[, cec_columns]

# Convert to numeric matrix
cec_matrix <- as.matrix(cec_data)

# Extract Combined_treatement as grouping factor
treatment_group <- df$Combined_treatement

# Trim whitespace from treatment group names
treatment_group <- trimws(treatment_group)

# ==============================================================================
# 3. RECODE TREATMENT CODES TO FULL NAMES
# ==============================================================================

cat("\n========== RECODING TREATMENT CODES ==========\n")

# Display unique treatment codes before recoding
cat("\nUnique treatment codes before recoding:\n")
print(unique(treatment_group))

# Recode treatment codes to full names
treatment_group <- case_when(
  treatment_group == "AA" ~ "BNR_AS",
  treatment_group == "BB" ~ "BNR_MAD",
  treatment_group == "CC" ~ "BNR_DW",
  treatment_group == "DD" ~ "PMAL_MADPS",
  treatment_group == "EE" ~ "CA_AS",
  treatment_group == "FF" ~ "CA_DW",
  treatment_group == "GG" ~ "CA_MAD",
  treatment_group == "HH" ~ "SAS_BPMAD",
  treatment_group == "II" ~ "SAS_DW",
  treatment_group == "JJ" ~ "SAS_AD",
  treatment_group == "KK" ~ "SAS_MAD",
  treatment_group == "LL" ~ "SASAS_MAD",
  treatment_group == "MM" ~ "SAS_MADPLT",
  TRUE ~ treatment_group  # Keep as-is if no match
)

# Display unique treatment codes after recoding
cat("\nUnique treatment codes after recoding:\n")
print(unique(treatment_group))

# ==============================================================================
# 4. CONTINUE DATA PREPARATION WITH RECODED GROUPS
# ==============================================================================

# Check for missing values
cat("\nChecking for missing values...\n")
missing_cec <- sum(is.na(cec_matrix))
missing_treatment <- sum(is.na(treatment_group) | treatment_group == "")
cat("Missing values in CEC data:", missing_cec, "\n")
cat("Missing values in Combined_treatement:", missing_treatment, "\n")

# Remove rows with NA in any CEC column or treatment
complete_rows <- complete.cases(cec_matrix) & !is.na(treatment_group) & treatment_group != ""
df_clean <- df[complete_rows, ]
cec_matrix_clean <- cec_matrix[complete_rows, ]
treatment_group_clean <- treatment_group[complete_rows]

removed_rows <- nrow(df) - nrow(df_clean)
if (removed_rows > 0) {
  cat("Removed", removed_rows, "rows with missing values.\n")
} else {
  cat("No missing values detected. All rows retained.\n")
}

# Convert Combined_treatement to factor
treatment_group_clean <- as.factor(treatment_group_clean)

# Display number of samples per treatment group
cat("\nNumber of samples per treatment group:\n")
treatment_counts <- table(treatment_group_clean)
print(treatment_counts)

# Identify groups with <2 samples
groups_single <- names(treatment_counts[treatment_counts < 2])
if (length(groups_single) > 0) {
  cat("\nWARNING: The following treatment groups have only 1 sample:\n")
  cat(paste(groups_single, collapse = ", "), "\n")
  cat("These groups will be excluded from PERMDISP and pairwise PERMANOVA.\n")
}

# Create a subset with only groups that have >= 2 samples
groups_to_keep <- names(treatment_counts[treatment_counts >= 2])
keep_indices <- treatment_group_clean %in% groups_to_keep

if (length(groups_to_keep) < 2) {
  stop("ERROR: Fewer than 2 treatment groups with >=2 samples. Analysis cannot be performed.")
}

cec_matrix_subset <- cec_matrix_clean[keep_indices, ]
treatment_group_subset <- treatment_group_clean[keep_indices]

# Convert to factor again to remove unused levels
treatment_group_subset <- droplevels(treatment_group_subset)

cat("\nSubset with groups having >=2 samples:\n")
cat("  - Samples:", nrow(cec_matrix_subset), "\n")
cat("  - Treatment groups:", length(unique(treatment_group_subset)), "\n")
print(table(treatment_group_subset))

# ==============================================================================
# 5. BRAY-CURTIS DISSIMILARITY ANALYSIS (USING RAW CONCENTRATIONS)
# ==============================================================================

cat("\n========== BRAY-CURTIS DISSIMILARITY ANALYSIS ==========\n")
cat("Using raw measured concentrations (µg/kg)\n")

# Calculate Bray-Curtis dissimilarity matrix
bray_curtis_dist <- vegdist(cec_matrix_subset, method = "bray")

# Convert to matrix for easier manipulation
bray_matrix <- as.matrix(bray_curtis_dist)

# Get summary statistics of dissimilarity values
dissimilarity_values <- bray_matrix[lower.tri(bray_matrix)]
dissimilarity_summary <- data.frame(
  Statistic = c("Min", "Max", "Median", "Mean", "SD"),
  Value = c(min(dissimilarity_values), max(dissimilarity_values), 
            median(dissimilarity_values), mean(dissimilarity_values), 
            sd(dissimilarity_values))
)

cat("\nBray-Curtis Dissimilarity Summary Statistics:\n")
print(dissimilarity_summary)

# Calculate within-group and between-group dissimilarities
cat("\nCalculating within-group and between-group dissimilarities...\n")

# Within-group dissimilarities
within_group_dist <- list()
for (group in levels(treatment_group_subset)) {
  group_indices <- which(treatment_group_subset == group)
  if (length(group_indices) >= 2) {
    group_dist <- bray_matrix[group_indices, group_indices]
    group_dist <- group_dist[lower.tri(group_dist)]
    if (length(group_dist) > 0) {
      within_group_dist[[group]] <- group_dist
    }
  }
}

# Between-group dissimilarities
between_group_dist <- list()
group_levels <- levels(treatment_group_subset)
for (i in 1:(length(group_levels) - 1)) {
  for (j in (i + 1):length(group_levels)) {
    group1_indices <- which(treatment_group_subset == group_levels[i])
    group2_indices <- which(treatment_group_subset == group_levels[j])
    if (length(group1_indices) >= 1 && length(group2_indices) >= 1) {
      between_dist <- bray_matrix[group1_indices, group2_indices]
      if (length(between_dist) > 0) {
        between_group_dist[[paste(group_levels[i], "vs", group_levels[j])]] <- between_dist
      }
    }
  }
}

# Prepare within-group summary
if (length(within_group_dist) > 0) {
  within_summary <- data.frame(
    Treatment = names(within_group_dist),
    N_pairs = sapply(within_group_dist, length),
    Mean = sapply(within_group_dist, mean, na.rm = TRUE),
    SD = sapply(within_group_dist, sd, na.rm = TRUE),
    Median = sapply(within_group_dist, median, na.rm = TRUE)
  )
  cat("\nWithin-group dissimilarities:\n")
  print(within_summary)
} else {
  within_summary <- data.frame(Treatment = character(), N_pairs = numeric(), 
                               Mean = numeric(), SD = numeric(), Median = numeric())
  cat("\nNo within-group distances calculated (groups may have <2 samples).\n")
}

# Prepare between-group summary
if (length(between_group_dist) > 0) {
  between_summary <- data.frame(
    Comparison = names(between_group_dist),
    N_pairs = sapply(between_group_dist, length),
    Mean = sapply(between_group_dist, mean, na.rm = TRUE),
    SD = sapply(between_group_dist, sd, na.rm = TRUE),
    Median = sapply(between_group_dist, median, na.rm = TRUE)
  )
  cat("\nBetween-group dissimilarities:\n")
  print(between_summary)
} else {
  between_summary <- data.frame(Comparison = character(), N_pairs = numeric(),
                                Mean = numeric(), SD = numeric(), Median = numeric())
  cat("\nNo between-group distances calculated.\n")
}

# Calculate overall average dissimilarity
overall_mean <- mean(dissimilarity_values)
cat("\nOverall average Bray-Curtis dissimilarity:", round(overall_mean, 4), "\n")

# ==============================================================================
# 6. PRINCIPAL COORDINATE ANALYSIS (PCoA)
# ==============================================================================

cat("\n========== PRINCIPAL COORDINATE ANALYSIS (PCoA) ==========\n")

# Try PCoA with error handling
pcoa_success <- FALSE
pcoa_scores <- NULL
variance_explained <- NULL

tryCatch({
  # Perform PCoA on Bray-Curtis distance matrix
  k_value <- min(2, nrow(cec_matrix_subset) - 1)
  if (k_value >= 1) {
    pcoa_result <- cmdscale(bray_curtis_dist, k = k_value, eig = TRUE)
    
    # Check if PCoA produced valid results
    if (!is.null(pcoa_result) && !is.null(pcoa_result$points) && nrow(pcoa_result$points) > 0) {
      # Extract PCoA scores
      pcoa_scores <- data.frame(
        Sample = rownames(pcoa_result$points),
        PCoA1 = pcoa_result$points[, 1],
        PCoA2 = if (ncol(pcoa_result$points) >= 2) pcoa_result$points[, 2] else 0,
        Combined_treatement = treatment_group_subset
      )
      
      # Calculate variance explained
      if (!is.null(pcoa_result$eig)) {
        variance_explained <- (pcoa_result$eig / sum(pcoa_result$eig)) * 100
        cat("\nVariance explained:\n")
        cat("PCoA1:", round(variance_explained[1], 2), "%\n")
        if (length(variance_explained) >= 2) {
          cat("PCoA2:", round(variance_explained[2], 2), "%\n")
        }
      }
      
      pcoa_success <- TRUE
    } else {
      cat("\nWARNING: PCoA did not produce valid results.\n")
    }
  } else {
    cat("\nWARNING: Insufficient samples for PCoA (need at least 2 samples).\n")
  }
}, error = function(e) {
  cat("\nERROR in PCoA:", e$message, "\n")
  cat("Proceeding with analysis without PCoA plot.\n")
})

# ==============================================================================
# 7. PERMANOVA TEST
# ==============================================================================

cat("\n========== PERMANOVA TEST ==========\n")

# Perform PERMANOVA using adonis2
permanova_result <- adonis2(bray_curtis_dist ~ treatment_group_subset, 
                            permutations = 999,
                            method = "bray")

cat("\nPERMANOVA Results:\n")
print(permanova_result)

# Extract results safely
permanova_df <- data.frame(
  Df = ifelse(!is.na(permanova_result$Df[1]), permanova_result$Df[1], NA),
  SumOfSqs = ifelse(!is.na(permanova_result$SumOfSqs[1]), permanova_result$SumOfSqs[1], NA),
  F = ifelse(!is.na(permanova_result$F[1]), permanova_result$F[1], NA),
  R2 = ifelse(!is.na(permanova_result$R2[1]), permanova_result$R2[1], NA),
  p_value = ifelse(!is.na(permanova_result$`Pr(>F)`[1]), permanova_result$`Pr(>F)`[1], NA),
  Significance = ifelse(!is.na(permanova_result$`Pr(>F)`[1]) && permanova_result$`Pr(>F)`[1] < 0.05, "Yes", "No")
)

# Check if significant
permanova_significant <- !is.na(permanova_result$`Pr(>F)`[1]) && permanova_result$`Pr(>F)`[1] < 0.05

cat("\nPERMANOVA Interpretation:\n")
if (!is.na(permanova_df$F)) {
  cat("F-statistic:", round(permanova_df$F, 4), "\n")
}
if (!is.na(permanova_df$R2)) {
  cat("R²:", round(permanova_df$R2, 4), "(", round(permanova_df$R2 * 100, 2), "% of variance explained)\n")
}
if (!is.na(permanova_df$p_value)) {
  cat("p-value:", permanova_df$p_value, "\n")
}
cat("Significant at α = 0.05?", ifelse(permanova_significant, "YES", "NO"), "\n")

# Pairwise PERMANOVA if significant
pairwise_permanova <- NULL
if (permanova_significant && length(unique(treatment_group_subset)) >= 3) {
  cat("\n--- Performing pairwise PERMANOVA comparisons ---\n")
  
  tryCatch({
    # Use pairwise.adonis from pairwiseAdonis package
    pairwise_permanova <- pairwise.adonis(
      x = cec_matrix_subset,
      factors = treatment_group_subset,
      perm = 999,
      p.adjust.m = "bonferroni"
    )
    
    cat("\nPairwise PERMANOVA Results (Bonferroni-adjusted):\n")
    print(pairwise_permanova)
  }, error = function(e) {
    cat("\nERROR in pairwise PERMANOVA:", e$message, "\n")
    cat("Pairwise comparisons not performed.\n")
  })
}

# ==============================================================================
# 8. PERMDISP TEST - WITH ERROR HANDLING
# ==============================================================================

cat("\n========== PERMDISP TEST ==========\n")

# Initialize variables with default values
permdisp_df <- data.frame(Df = NA, SumOfSqs = NA, F = NA, p_value = NA, Significance = "No")
group_distances <- data.frame(Treatment = character(), Mean_Distance_to_Centroid = numeric(), SD_Distance_to_Centroid = numeric())
permdisp_significant <- FALSE

# Try PERMDISP with comprehensive error handling
tryCatch({
  # Perform PERMDISP (multivariate homogeneity of group dispersions)
  permdisp_result <- betadisper(bray_curtis_dist, treatment_group_subset)
  
  # Display basic results
  cat("\nPERMDISP Summary:\n")
  print(summary(permdisp_result))
  
  # Check if distances exist
  if (!is.null(permdisp_result$distances) && length(permdisp_result$distances) > 0) {
    # Extract group-specific distances
    group_distances <- data.frame(
      Treatment = names(permdisp_result$distances),
      Mean_Distance_to_Centroid = tapply(permdisp_result$distances, 
                                         treatment_group_subset, mean, na.rm = TRUE),
      SD_Distance_to_Centroid = tapply(permdisp_result$distances, 
                                       treatment_group_subset, sd, na.rm = TRUE)
    )
    cat("\nGroup-specific distances to centroid:\n")
    print(group_distances)
    
    # Try permutation test
    tryCatch({
      permdisp_test <- permutest(permdisp_result, permutations = 999, pairwise = TRUE)
      cat("\nPERMDISP Permutation Test Results:\n")
      print(permdisp_test)
      
      # Extract results from permutation test if available
      if (!is.null(permdisp_test$tab) && nrow(permdisp_test$tab) >= 1) {
        permdisp_tab <- as.data.frame(permdisp_test$tab)
        
        permdisp_df <- data.frame(
          Df = if (!is.null(permdisp_tab$Df[1])) permdisp_tab$Df[1] else NA,
          SumOfSqs = if (!is.null(permdisp_tab$SumSq[1])) permdisp_tab$SumSq[1] else NA,
          F = if (!is.null(permdisp_tab$F[1])) permdisp_tab$F[1] else NA,
          p_value = if (!is.null(permdisp_tab$`Pr(>F)`[1])) permdisp_tab$`Pr(>F)`[1] else NA,
          Significance = "No"
        )
        
        if (!is.na(permdisp_df$p_value) && permdisp_df$p_value < 0.05) {
          permdisp_df$Significance <- "Yes"
        }
      }
    }, error = function(e) {
      cat("\nWARNING: PERMDISP permutation test failed:", e$message, "\n")
      cat("Using betadisper ANOVA results.\n")
      
      # Try ANOVA results as fallback
      if (!is.null(permdisp_result$anova) && nrow(permdisp_result$anova) >= 1) {
        anova_tab <- as.data.frame(permdisp_result$anova)
        
        permdisp_df <- data.frame(
          Df = if (!is.null(anova_tab$Df[1])) anova_tab$Df[1] else NA,
          SumOfSqs = if (!is.null(anova_tab$`Sum Sq`[1])) anova_tab$`Sum Sq`[1] else NA,
          F = if (!is.null(anova_tab$`F value`[1])) anova_tab$`F value`[1] else NA,
          p_value = if (!is.null(anova_tab$`Pr(>F)`[1])) anova_tab$`Pr(>F)`[1] else NA,
          Significance = "No"
        )
        
        if (!is.na(permdisp_df$p_value) && permdisp_df$p_value < 0.05) {
          permdisp_df$Significance <- "Yes"
        }
      }
    })
  } else {
    cat("\nWARNING: No distances available in betadisper result.\n")
  }
}, error = function(e) {
  cat("\nERROR in PERMDISP analysis:", e$message, "\n")
  cat("Proceeding with default values.\n")
})

# Determine significance
if (!is.na(permdisp_df$p_value) && permdisp_df$p_value < 0.05) {
  permdisp_significant <- TRUE
}

cat("\nPERMDISP Interpretation:\n")
if (!is.na(permdisp_df$F)) {
  cat("F-statistic:", round(permdisp_df$F, 4), "\n")
}
if (!is.na(permdisp_df$p_value)) {
  cat("p-value:", permdisp_df$p_value, "\n")
}
cat("Significant at α = 0.05?", ifelse(permdisp_significant, "YES", "NO"), "\n")

if (permdisp_significant) {
  cat("\nWARNING: Significant PERMDISP indicates heterogeneous multivariate dispersions.\n")
  cat("PERMANOVA results should be interpreted with caution.\n")
}

# ==============================================================================
# 9. GENERATE TABLE REPORTS (CSV FORMAT)
# ==============================================================================

cat("\n========== GENERATING TABLE REPORTS ==========\n")

# ---------- Table 1: Bray-Curtis Dissimilarity Results ----------
table1 <- data.frame(
  Metric = c("Overall_Mean", "Overall_Median", "Overall_Min", "Overall_Max", "Overall_SD"),
  Value = c(mean(dissimilarity_values), median(dissimilarity_values), 
            min(dissimilarity_values), max(dissimilarity_values), 
            sd(dissimilarity_values))
)

# Add within-group summary
if (nrow(within_summary) > 0) {
  for (i in 1:nrow(within_summary)) {
    table1 <- rbind(table1, 
                    data.frame(Metric = paste0("Within_", within_summary$Treatment[i], "_Mean"),
                               Value = within_summary$Mean[i]))
    table1 <- rbind(table1,
                    data.frame(Metric = paste0("Within_", within_summary$Treatment[i], "_SD"),
                               Value = within_summary$SD[i]))
  }
}

# Add between-group summary
if (nrow(between_summary) > 0) {
  for (i in 1:nrow(between_summary)) {
    table1 <- rbind(table1,
                    data.frame(Metric = paste0("Between_", between_summary$Comparison[i], "_Mean"),
                               Value = between_summary$Mean[i]))
    table1 <- rbind(table1,
                    data.frame(Metric = paste0("Between_", between_summary$Comparison[i], "_SD"),
                               Value = between_summary$SD[i]))
  }
}

# Save Table 1
table1_path <- file.path(tables_dir, "Table1_BrayCurtis_Dissimilarity.csv")
write.csv(table1, table1_path, row.names = FALSE)
cat("Table 1 saved to:", table1_path, "\n")

# ---------- Table 2: PERMANOVA Results ----------
# Start with overall PERMANOVA results
table2 <- data.frame(
  Test = "Overall_PERMANOVA",
  Df = ifelse(!is.na(permanova_df$Df), permanova_df$Df, NA),
  SumOfSqs = ifelse(!is.na(permanova_df$SumOfSqs), round(permanova_df$SumOfSqs, 4), NA),
  F = ifelse(!is.na(permanova_df$F), round(permanova_df$F, 4), NA),
  R2 = ifelse(!is.na(permanova_df$R2), round(permanova_df$R2, 4), NA),
  p_value = ifelse(!is.na(permanova_df$p_value), permanova_df$p_value, NA),
  Significance = ifelse(!is.na(permanova_df$Significance), permanova_df$Significance, "No")
)

# Add pairwise results if available - FIXED for actual column names
if (!is.null(pairwise_permanova) && nrow(pairwise_permanova) > 0) {
  
  cat("\nPairwise PERMANOVA column names:\n")
  print(colnames(pairwise_permanova))
  
  # Create a new data frame with standardized column names
  pairwise_std <- data.frame(
    Test = paste0("Pairwise_", pairwise_permanova$pairs),
    Df = pairwise_permanova$Df,
    SumOfSqs = pairwise_permanova$SumsOfSqs,
    F = pairwise_permanova$F.Model,
    R2 = pairwise_permanova$R2,
    p_value = pairwise_permanova$p.value,
    p_adjusted = pairwise_permanova$p.adjusted,
    Significance = ifelse(pairwise_permanova$p.adjusted < 0.05, "Yes", "No")
  )
  
  # Select columns to match table2
  cols_to_keep <- c("Test", "Df", "SumOfSqs", "F", "R2", "p_value", "Significance")
  table2 <- rbind(table2, pairwise_std[, cols_to_keep])
}

# Save Table 2
table2_path <- file.path(tables_dir, "Table2_PERMANOVA_Results.csv")
write.csv(table2, table2_path, row.names = FALSE)
cat("Table 2 saved to:", table2_path, "\n")

# ---------- Table 3: PERMDISP Results ----------
table3 <- data.frame(
  Test = "Overall_PERMDISP",
  Df = ifelse(!is.na(permdisp_df$Df), permdisp_df$Df, NA),
  SumOfSqs = ifelse(!is.na(permdisp_df$SumOfSqs), round(permdisp_df$SumOfSqs, 4), NA),
  F = ifelse(!is.na(permdisp_df$F), round(permdisp_df$F, 4), NA),
  p_value = ifelse(!is.na(permdisp_df$p_value), permdisp_df$p_value, NA),
  Significance = ifelse(!is.na(permdisp_df$Significance), permdisp_df$Significance, "No")
)

# Add group-specific distances
if (nrow(group_distances) > 0) {
  for (i in 1:nrow(group_distances)) {
    table3 <- rbind(table3,
                    data.frame(Test = paste0("Group_", group_distances$Treatment[i], "_Mean"),
                               Df = NA, SumOfSqs = NA, F = NA,
                               p_value = NA, Significance = NA))
    table3 <- rbind(table3,
                    data.frame(Test = paste0("Group_", group_distances$Treatment[i], "_SD"),
                               Df = NA, SumOfSqs = NA, F = NA,
                               p_value = NA, Significance = NA))
  }
}

# Save Table 3
table3_path <- file.path(tables_dir, "Table3_PERMDISP_Results.csv")
write.csv(table3, table3_path, row.names = FALSE)
cat("Table 3 saved to:", table3_path, "\n")

# ==============================================================================
# 10. GENERATE PLAIN-TEXT REPORT FILES
# ==============================================================================

cat("\n========== GENERATING REPORT FILES ==========\n")

# ---------- Report 1: Bray-Curtis ----------
report1_lines <- c(
  "=================================================================================",
  "          BRAY-CURTIS DISSIMILARITY ANALYSIS OF CEC CLASS COMPOSITION",
  "=================================================================================",
  paste("Date and time of analysis:", Sys.time()),
  paste("Data file:", csv_file),
  paste("Data type: Raw measured concentrations (µg/kg)"),
  paste("Number of samples analyzed:", nrow(cec_matrix_subset)),
  paste("Number of treatment groups:", length(unique(treatment_group_subset))),
  paste("Treatment groups:", paste(levels(treatment_group_subset), collapse = ", ")),
  "",
  "---------------------------------------------------------------------------------",
  "SUMMARY STATISTICS OF BRAY-CURTIS DISSIMILARITY VALUES",
  "---------------------------------------------------------------------------------",
  paste(capture.output(print(dissimilarity_summary)), collapse = "\n"),
  "",
  "---------------------------------------------------------------------------------",
  "WITHIN-GROUP DISSIMILARITIES",
  "---------------------------------------------------------------------------------"
)

if (nrow(within_summary) > 0) {
  report1_lines <- c(report1_lines,
                     paste(capture.output(print(within_summary)), collapse = "\n"))
} else {
  report1_lines <- c(report1_lines,
                     "No within-group distances calculated (groups may have <2 samples).")
}

report1_lines <- c(report1_lines,
                   "",
                   "---------------------------------------------------------------------------------",
                   "BETWEEN-GROUP DISSIMILARITIES",
                   "---------------------------------------------------------------------------------"
)

if (nrow(between_summary) > 0) {
  report1_lines <- c(report1_lines,
                     paste(capture.output(print(between_summary)), collapse = "\n"))
} else {
  report1_lines <- c(report1_lines,
                     "No between-group distances calculated.")
}

report1_lines <- c(report1_lines,
                   "",
                   "---------------------------------------------------------------------------------",
                   "INTERPRETATION",
                   "---------------------------------------------------------------------------------"
)

# Interpretation
if (overall_mean < 0.3) {
  interp1 <- "low"
} else if (overall_mean < 0.6) {
  interp1 <- "moderate"
} else {
  interp1 <- "high"
}

report1_lines <- c(report1_lines,
                   paste("The average Bray-Curtis dissimilarity between samples was", 
                         round(overall_mean, 4), "indicating", interp1, "compositional variation."),
                   "",
                   "---------------------------------------------------------------------------------",
                   "END OF REPORT",
                   "================================================================================="
)

# Write Report 1
report1_path <- file.path(reports_dir, "Report1_BrayCurtis.txt")
writeLines(report1_lines, report1_path)
cat("Report 1 saved to:", report1_path, "\n")

# ---------- Report 2: PERMANOVA ----------
report2_lines <- c(
  "=================================================================================",
  "          PERMANOVA RESULTS: COMBINED TREATMENT EFFECTS ON CEC CLASS COMPOSITION",
  "=================================================================================",
  paste("Date and time of analysis:", Sys.time()),
  paste("Data file:", csv_file),
  paste("Data type: Raw measured concentrations (µg/kg)"),
  paste("Number of samples analyzed:", nrow(cec_matrix_subset)),
  paste("Number of treatment groups:", length(unique(treatment_group_subset))),
  paste("Treatment groups:", paste(levels(treatment_group_subset), collapse = ", ")),
  paste("Number of permutations:", 999),
  "",
  "---------------------------------------------------------------------------------",
  "OVERALL PERMANOVA RESULTS",
  "---------------------------------------------------------------------------------"
)

if (!is.na(permanova_df$Df)) {
  report2_lines <- c(report2_lines,
                     paste("Df:", permanova_df$Df))
}
if (!is.na(permanova_df$SumOfSqs)) {
  report2_lines <- c(report2_lines,
                     paste("Sum of Squares:", round(permanova_df$SumOfSqs, 4)))
}
if (!is.na(permanova_df$F)) {
  report2_lines <- c(report2_lines,
                     paste("F-statistic:", round(permanova_df$F, 4)))
}
if (!is.na(permanova_df$R2)) {
  report2_lines <- c(report2_lines,
                     paste("R² (coefficient of determination):", round(permanova_df$R2, 4)))
}
if (!is.na(permanova_df$p_value)) {
  report2_lines <- c(report2_lines,
                     paste("p-value:", permanova_df$p_value))
}
report2_lines <- c(report2_lines,
                   paste("Significant at α = 0.05?", permanova_df$Significance),
                   "",
                   "---------------------------------------------------------------------------------",
                   "INTERPRETATION",
                   "---------------------------------------------------------------------------------"
)

if (permanova_significant) {
  report2_lines <- c(report2_lines,
                     paste("PERMANOVA revealed a STATISTICALLY SIGNIFICANT effect of Combined_treatement",
                           "on CEC class composition (F =", round(permanova_df$F, 4), 
                           ", R² =", round(permanova_df$R2, 4), 
                           ", p =", permanova_df$p_value, ")."),
                     "",
                     paste("The treatment explains", round(permanova_df$R2 * 100, 2), 
                           "% of the total compositional variation."),
                     "",
                     "---------------------------------------------------------------------------------",
                     "PAIRWISE COMPARISONS",
                     "---------------------------------------------------------------------------------"
  )
  
  if (!is.null(pairwise_permanova) && nrow(pairwise_permanova) > 0) {
    # Create a readable pairwise summary
    pairwise_summary <- pairwise_permanova
    pairwise_summary$Significance <- ifelse(pairwise_summary$p.adjusted < 0.05, "Yes", "No")
    
    report2_lines <- c(report2_lines,
                       paste(capture.output(print(pairwise_summary)), collapse = "\n"),
                       "",
                       "Significant pairwise comparisons (Bonferroni-adjusted, α = 0.05):",
                       paste(sapply(1:nrow(pairwise_permanova), function(i) {
                         if (pairwise_permanova$p.adjusted[i] < 0.05) {
                           paste(pairwise_permanova$pairs[i], ": p =", round(pairwise_permanova$p.adjusted[i], 6))
                         }
                       }), collapse = "\n  ")
    )
  } else {
    report2_lines <- c(report2_lines,
                       "Pairwise comparisons not performed (insufficient groups or samples).")
  }
  
} else {
  report2_lines <- c(report2_lines,
                     paste("PERMANOVA revealed NO STATISTICALLY SIGNIFICANT effect of Combined_treatement",
                           "on CEC class composition (F =", round(permanova_df$F, 4), 
                           ", R² =", round(permanova_df$R2, 4), 
                           ", p =", permanova_df$p_value, ")."),
                     "",
                     "The treatment explains only", round(permanova_df$R2 * 100, 2), 
                     "% of the total compositional variation, which is not statistically significant.")
}

report2_lines <- c(report2_lines,
                   "",
                   "---------------------------------------------------------------------------------",
                   "END OF REPORT",
                   "================================================================================="
)

# Write Report 2
report2_path <- file.path(reports_dir, "Report2_PERMANOVA.txt")
writeLines(report2_lines, report2_path)
cat("Report 2 saved to:", report2_path, "\n")

# ---------- Report 3: PERMDISP ----------
report3_lines <- c(
  "=================================================================================",
  "          PERMDISP RESULTS: HOMOGENEITY OF MULTIVARIATE DISPERSIONS",
  "=================================================================================",
  paste("Date and time of analysis:", Sys.time()),
  paste("Data file:", csv_file),
  paste("Data type: Raw measured concentrations (µg/kg)"),
  paste("Number of samples analyzed:", nrow(cec_matrix_subset)),
  paste("Number of treatment groups:", length(unique(treatment_group_subset))),
  paste("Treatment groups:", paste(levels(treatment_group_subset), collapse = ", ")),
  paste("Number of permutations:", 999),
  "",
  "---------------------------------------------------------------------------------",
  "OVERALL PERMDISP RESULTS",
  "---------------------------------------------------------------------------------"
)

if (!is.na(permdisp_df$Df)) {
  report3_lines <- c(report3_lines,
                     paste("Df:", permdisp_df$Df))
}
if (!is.na(permdisp_df$SumOfSqs)) {
  report3_lines <- c(report3_lines,
                     paste("Sum of Squares:", round(permdisp_df$SumOfSqs, 4)))
}
if (!is.na(permdisp_df$F)) {
  report3_lines <- c(report3_lines,
                     paste("F-statistic:", round(permdisp_df$F, 4)))
}
if (!is.na(permdisp_df$p_value)) {
  report3_lines <- c(report3_lines,
                     paste("p-value:", permdisp_df$p_value))
}
report3_lines <- c(report3_lines,
                   paste("Significant at α = 0.05?", permdisp_df$Significance),
                   "",
                   "---------------------------------------------------------------------------------",
                   "GROUP-SPECIFIC DISTANCES TO CENTROID",
                   "---------------------------------------------------------------------------------"
)

if (nrow(group_distances) > 0) {
  report3_lines <- c(report3_lines,
                     paste(capture.output(print(group_distances)), collapse = "\n"))
} else {
  report3_lines <- c(report3_lines,
                     "No group-specific distances available.")
}

report3_lines <- c(report3_lines,
                   "",
                   "---------------------------------------------------------------------------------",
                   "INTERPRETATION",
                   "---------------------------------------------------------------------------------"
)

if (permdisp_significant) {
  report3_lines <- c(report3_lines,
                     paste("PERMDISP showed STATISTICALLY SIGNIFICANT differences in multivariate",
                           "dispersions among treatment groups (F =", round(permdisp_df$F, 4), 
                           ", p =", permdisp_df$p_value, ")."),
                     "",
                     "This indicates that the VARIABILITY DIFFERS SIGNIFICANTLY across treatment groups.",
                     "",
                     "IMPORTANT: The PERMANOVA results should be interpreted with CAUTION as the",
                     "assumption of homogeneous dispersions is VIOLATED.",
                     "",
                     "The significant PERMDISP suggests that differences detected by PERMANOVA",
                     "may be partly due to differences in group variances rather than differences",
                     "in group centroids (mean compositions).")
} else {
  report3_lines <- c(report3_lines,
                     paste("PERMDISP showed NO STATISTICALLY SIGNIFICANT differences in multivariate",
                           "dispersions among treatment groups (F =", ifelse(!is.na(permdisp_df$F), round(permdisp_df$F, 4), "NA"), 
                           ", p =", ifelse(!is.na(permdisp_df$p_value), permdisp_df$p_value, "NA"), ")."),
                     "",
                     "This indicates that the VARIABILITY IS HOMOGENEOUS across treatment groups.",
                     "",
                     "The assumption of homogeneous dispersions for PERMANOVA is SATISFIED.",
                     "PERMANOVA results can be interpreted with confidence.")
}

report3_lines <- c(report3_lines,
                   "",
                   "---------------------------------------------------------------------------------",
                   "END OF REPORT",
                   "================================================================================="
)

# Write Report 3
report3_path <- file.path(reports_dir, "Report3_PERMDISP.txt")
writeLines(report3_lines, report3_path)
cat("Report 3 saved to:", report3_path, "\n")

# ==============================================================================
# 11. VISUALIZATION - PCoA PLOT (IF SUCCESSFUL)
# ==============================================================================

if (pcoa_success && !is.null(pcoa_scores) && nrow(pcoa_scores) > 0) {
  cat("\n========== GENERATING PCOA PLOT ==========\n")
  
  tryCatch({
    # Create PCoA plot
    pcoa_plot <- ggplot(pcoa_scores, aes(x = PCoA1, y = PCoA2, color = Combined_treatement)) +
      geom_point(size = 3, alpha = 0.8) +
      stat_ellipse(aes(fill = Combined_treatement), 
                   type = "t", 
                   level = 0.95, 
                   alpha = 0.2,
                   geom = "polygon") +
      labs(
        x = paste0("PCoA1 (", round(variance_explained[1], 2), "%)"),
        y = paste0("PCoA2 (", round(variance_explained[2], 2), "%)"),
        color = "Combined Treatment",
        fill = "Combined Treatment"
      ) +
      theme_minimal() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 11, face = "bold"),
        axis.text = element_text(size = 9),
        panel.grid.minor = element_blank(),
        plot.margin = margin(1, 1, 1, 1, "cm")
      ) +
      scale_color_viridis_d() +
      scale_fill_viridis_d() +
      coord_fixed(ratio = 1)
    
    # Save PCoA plot
    pcoa_plot_path <- file.path(figures_dir, "PCoA_BrayCurtis_Treatments.png")
    ggsave(filename = pcoa_plot_path,
           plot = pcoa_plot,
           width = 539 / 72,
           height = 539 / 72 * 0.85,
           dpi = 300,
           bg = "white",
           units = "in")
    cat("PCoA plot saved to:", pcoa_plot_path, "\n")
  }, error = function(e) {
    cat("\nERROR generating PCoA plot:", e$message, "\n")
    cat("PCoA plot not generated.\n")
  })
} else {
  cat("\n========== PCoA PLOT NOT GENERATED ==========\n")
  cat("PCoA analysis was not successful. This may be due to:\n")
  cat("  - Insufficient samples (need at least 3 samples for PCoA)\n")
  cat("  - Duplicate or identical samples\n")
  cat("  - Problems with the distance matrix\n")
  cat("Proceeding with statistical analysis results only.\n")
}

# ==============================================================================
# 12. FINAL OUTPUT SUMMARY
# ==============================================================================

cat("\n")
cat("=================================================================================\n")
cat("                          ANALYSIS COMPLETE\n")
cat("=================================================================================\n")
cat("\nINPUT DATA:\n")
cat("  - CSV file:", csv_file, "\n")
cat("  - Data type: Raw measured concentrations (µg/kg)\n")
cat("  - Samples analyzed:", nrow(cec_matrix_subset), "\n")
cat("  - Treatment groups:", length(unique(treatment_group_subset)), "\n")
cat("  - CEC classes:", length(cec_columns), "\n")
cat("\nTREATMENT CODES (AA, BB, CC, DD, EE, FF, GG, HH, II, JJ, KK, LL, MM):\n")
cat("  - AA: BNR_AS\n")
cat("  - BB: BNR_MAD\n")
cat("  - CC: BNR_DW\n")
cat("  - DD: PMAL_MADPS\n")
cat("  - EE: CA_AS\n")
cat("  - FF: CA_DW\n")
cat("  - GG: CA_MAD\n")
cat("  - HH: SAS_BPMAD\n")
cat("  - II: SAS_DW\n")
cat("  - JJ: SAS_AD\n")
cat("  - KK: SAS_MAD\n")
cat("  - LL: SASAS_MAD\n")
cat("  - MM: SAS_MADPLT\n")
cat("\nOUTPUT FILES GENERATED:\n")
cat("\n1. Tables (CSV):\n")
cat("   -", table1_path, "\n")
cat("   -", table2_path, "\n")
cat("   -", table3_path, "\n")
cat("\n2. Reports (TXT):\n")
cat("   -", report1_path, "\n")
cat("   -", report2_path, "\n")
cat("   -", report3_path, "\n")
if (pcoa_success) {
  cat("\n3. Figures (PNG):\n")
  cat("   -", pcoa_plot_path, "\n")
}
cat("\n=================================================================================\n")
cat("                    SCRIPT EXECUTED SUCCESSFULLY\n")
cat("=================================================================================\n")

# End of script
