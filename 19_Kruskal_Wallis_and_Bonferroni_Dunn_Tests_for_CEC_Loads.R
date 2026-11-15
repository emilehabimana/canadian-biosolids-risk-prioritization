################################################################################
# SCRIPT: Kruskal-Wallis and Bonferroni-Dunn Tests for CEC Loads
# Objective: Compare CEC loads across WWTP size and Community size categories
# Category codes: SM=Small, MED=Medium, LRG=Large, VL=Very Large
#                S=Small, M=Medium, L=Large, MET=Metropolitan
# NOTE: Script automatically trims whitespace from category codes
################################################################################

# ==============================================================================
# 0. SETUP - Load Libraries and Configure Environment
# ==============================================================================

# Load required packages
library(dplyr)
library(rstatix)
library(FSA)
library(tidyr)

# Disable scientific notation for cleaner output
options(scipen = 999)

# Set seed for reproducibility
set.seed(123)

# ==============================================================================
# 1. DATA IMPORT AND PREPARATION
# ==============================================================================

# Define file path
file_path <- "C:/Users/Admin/Documents/Data/loads_and_fluxes.csv"

# Check if file exists
if (!file.exists(file_path)) {
  stop(paste("ERROR: File not found at:", file_path, 
             "\nPlease verify the path and file name."))
}

# Read CSV file
df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

# Display structure and first few rows for verification
cat("\n========== DATA IMPORT VERIFICATION ==========\n")
cat("File successfully loaded from:", file_path, "\n")
cat("\nData frame dimensions:", dim(df)[1], "rows x", dim(df)[2], "columns\n")
cat("\nColumn names:\n")
print(colnames(df))
cat("\nFirst 6 rows of data:\n")
print(head(df))

# ==============================================================================
# 2. TRIM WHITESPACE AND RECODE CATEGORICAL VARIABLES
# ==============================================================================

cat("\n========== TRIMMING WHITESPACE AND RECODING ==========\n")

# IMPORTANT: Trim whitespace from categorical columns FIRST
df$WWTP_size <- trimws(df$WWTP_size)
df$Community_size <- trimws(df$Community_size)

# Display original values after trimming (should now match expected codes)
cat("\nWWTP_size values after trimming whitespace:\n")
print(unique(df$WWTP_size))
cat("\nCommunity_size values after trimming whitespace:\n")
print(unique(df$Community_size))

# Recode WWTP_size from codes to full names
df$WWTP_size <- case_when(
  df$WWTP_size == "SM" ~ "Small",
  df$WWTP_size == "MED" ~ "Medium",
  df$WWTP_size == "LRG" ~ "Large",
  df$WWTP_size == "VL" ~ "Very Large",
  TRUE ~ as.character(df$WWTP_size)  # Keep as-is if no match
)

# Recode Community_size from codes to full names
df$Community_size <- case_when(
  df$Community_size == "S" ~ "Small",
  df$Community_size == "M" ~ "Medium",
  df$Community_size == "L" ~ "Large",
  df$Community_size == "MET" ~ "Metropolitan",
  TRUE ~ as.character(df$Community_size)  # Keep as-is if no match
)

# Check for any values that didn't match the expected codes (after trimming)
wwtp_unmatched <- unique(df$WWTP_size[!df$WWTP_size %in% c("Small", "Medium", "Large", "Very Large")])
community_unmatched <- unique(df$Community_size[!df$Community_size %in% c("Small", "Medium", "Large", "Metropolitan")])

if (length(wwtp_unmatched) > 0) {
  cat("\nWARNING: These WWTP_size values did not match expected codes (SM, MED, LRG, VL):\n")
  print(wwtp_unmatched)
  cat("Please check for typos or extra characters.\n")
} else {
  cat("\nAll WWTP_size values successfully matched expected codes.\n")
}

if (length(community_unmatched) > 0) {
  cat("\nWARNING: These Community_size values did not match expected codes (S, M, L, MET):\n")
  print(community_unmatched)
  cat("Please check for typos or extra characters.\n")
} else {
  cat("All Community_size values successfully matched expected codes.\n")
}

# Convert to ordered factors with proper levels
df$WWTP_size <- factor(df$WWTP_size, 
                       levels = c("Small", "Medium", "Large", "Very Large"),
                       ordered = TRUE)

df$Community_size <- factor(df$Community_size,
                            levels = c("Small", "Medium", "Large", "Metropolitan"),
                            ordered = TRUE)

# Verify the conversion
cat("\n========== CATEGORY CONVERSION VERIFICATION ==========\n")
cat("\nWWTP_size after recoding:\n")
print(table(df$WWTP_size))
cat("\nCommunity_size after recoding:\n")
print(table(df$Community_size))

# ==============================================================================
# 3. CHECK FOR MISSING VALUES
# ==============================================================================

cat("\n========== MISSING DATA CHECK ==========\n")
missing_loads <- sum(is.na(df$CEC_loads))
missing_wwtp <- sum(is.na(df$WWTP_size) | df$WWTP_size == "")
missing_community <- sum(is.na(df$Community_size) | df$Community_size == "")
cat("Missing values in CEC_loads:", missing_loads, "\n")
cat("Missing values in WWTP_size:", missing_wwtp, "\n")
cat("Missing values in Community_size:", missing_community, "\n")

# Remove rows with NA in critical columns
df_clean <- df %>%
  filter(!is.na(CEC_loads) & 
           !is.na(WWTP_size) & WWTP_size != "" &
           !is.na(Community_size) & Community_size != "")

removed_rows <- nrow(df) - nrow(df_clean)
if (removed_rows > 0) {
  cat("Removed", removed_rows, "rows with missing values.\n")
} else {
  cat("No missing values detected. All rows retained.\n")
}

# Check if we have at least 2 groups for each variable (required for Kruskal-Wallis)
wwtp_groups <- unique(df_clean$WWTP_size)
community_groups <- unique(df_clean$Community_size)

cat("\nNumber of WWTP_size groups:", length(wwtp_groups), "\n")
cat("WWTP_size groups:", paste(wwtp_groups, collapse = ", "), "\n")
cat("Number of Community_size groups:", length(community_groups), "\n")
cat("Community_size groups:", paste(community_groups, collapse = ", "), "\n")

if (length(wwtp_groups) < 2) {
  stop("ERROR: WWTP_size has fewer than 2 groups. Kruskal-Wallis test cannot be performed.")
}

if (length(community_groups) < 2) {
  stop("ERROR: Community_size has fewer than 2 groups. Kruskal-Wallis test cannot be performed.")
}

# ==============================================================================
# 4. GENERATE SUMMARY STATISTICS
# ==============================================================================

# Summary statistics grouped by WWTP_size
cat("\n========== SUMMARY STATISTICS BY WWTP SIZE ==========\n")
summary_wwtp <- df_clean %>%
  group_by(WWTP_size) %>%
  summarise(
    n = n(),
    Mean = mean(CEC_loads, na.rm = TRUE),
    SD = sd(CEC_loads, na.rm = TRUE),
    Median = median(CEC_loads, na.rm = TRUE),
    IQR = IQR(CEC_loads, na.rm = TRUE),
    Min = min(CEC_loads, na.rm = TRUE),
    Max = max(CEC_loads, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_wwtp)

# Summary statistics grouped by Community_size
cat("\n========== SUMMARY STATISTICS BY COMMUNITY SIZE ==========\n")
summary_community <- df_clean %>%
  group_by(Community_size) %>%
  summarise(
    n = n(),
    Mean = mean(CEC_loads, na.rm = TRUE),
    SD = sd(CEC_loads, na.rm = TRUE),
    Median = median(CEC_loads, na.rm = TRUE),
    IQR = IQR(CEC_loads, na.rm = TRUE),
    Min = min(CEC_loads, na.rm = TRUE),
    Max = max(CEC_loads, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_community)

# Combined summary: WWTP_size × Community_size interaction
cat("\n========== COMBINED SUMMARY (WWTP_SIZE × COMMUNITY_SIZE) ==========\n")
combined_summary <- df_clean %>%
  group_by(WWTP_size, Community_size) %>%
  summarise(
    n = n(),
    Median = median(CEC_loads, na.rm = TRUE),
    IQR = IQR(CEC_loads, na.rm = TRUE),
    .groups = "drop"
  )
print(combined_summary)

# ==============================================================================
# 5. FIRST KRUSKAL-WALLIS TEST: CEC_loads vs. WWTP_size
# ==============================================================================

cat("\n========== TEST 1: KRUSKAL-WALLIS - CEC LOADS vs. WWTP SIZE ==========\n")

# Perform Kruskal-Wallis test
kw_test1 <- kruskal.test(CEC_loads ~ WWTP_size, data = df_clean)

# Extract results
chi_sq1 <- kw_test1$statistic
df1 <- kw_test1$parameter
p_val1 <- kw_test1$p.value
significant1 <- p_val1 < 0.05

# Display results
cat("Chi-squared:", chi_sq1, "\n")
cat("Degrees of freedom:", df1, "\n")
cat("p-value:", p_val1, "\n")
cat("Significant at α = 0.05?", ifelse(significant1, "YES", "NO"), "\n")

# Post-hoc test if significant
posthoc1 <- NULL
if (significant1) {
  cat("\n--- Significant difference detected. Performing Bonferroni-Dunn post-hoc test ---\n")
  
  # Use dunnTest from FSA package
  posthoc1 <- dunnTest(CEC_loads ~ WWTP_size, 
                       data = df_clean, 
                       method = "bonferroni")
  
  # Extract and format post-hoc results
  ph_results1 <- posthoc1$res
  ph_results1$Significant <- ph_results1$P.adj < 0.05
  
  cat("\nPost-hoc pairwise comparisons (Bonferroni-adjusted):\n")
  print(ph_results1)
  
} else {
  cat("\n--- No significant difference detected. Post-hoc test not performed. ---\n")
}

# ==============================================================================
# 6. SECOND KRUSKAL-WALLIS TEST: CEC_loads vs. Community_size
# ==============================================================================

cat("\n========== TEST 2: KRUSKAL-WALLIS - CEC LOADS vs. COMMUNITY SIZE ==========\n")

# Perform Kruskal-Wallis test
kw_test2 <- kruskal.test(CEC_loads ~ Community_size, data = df_clean)

# Extract results
chi_sq2 <- kw_test2$statistic
df2 <- kw_test2$parameter
p_val2 <- kw_test2$p.value
significant2 <- p_val2 < 0.05

# Display results
cat("Chi-squared:", chi_sq2, "\n")
cat("Degrees of freedom:", df2, "\n")
cat("p-value:", p_val2, "\n")
cat("Significant at α = 0.05?", ifelse(significant2, "YES", "NO"), "\n")

# Post-hoc test if significant
posthoc2 <- NULL
if (significant2) {
  cat("\n--- Significant difference detected. Performing Bonferroni-Dunn post-hoc test ---\n")
  
  # Use dunnTest from FSA package
  posthoc2 <- dunnTest(CEC_loads ~ Community_size, 
                       data = df_clean, 
                       method = "bonferroni")
  
  # Extract and format post-hoc results
  ph_results2 <- posthoc2$res
  ph_results2$Significant <- ph_results2$P.adj < 0.05
  
  cat("\nPost-hoc pairwise comparisons (Bonferroni-adjusted):\n")
  print(ph_results2)
  
} else {
  cat("\n--- No significant difference detected. Post-hoc test not performed. ---\n")
}

# ==============================================================================
# 7. GENERATE TABLE REPORTS (CSV FORMAT)
# ==============================================================================

# Define output directories
tables_dir <- "C:/Users/Admin/Documents/Article_2/Data/Tables"
reports_dir <- "C:/Users/Admin/Documents/Article_2/Data/Reports"

# Create directories if they don't exist
if (!dir.exists(tables_dir)) {
  dir.create(tables_dir, recursive = TRUE)
  cat("\nCreated directory:", tables_dir, "\n")
}
if (!dir.exists(reports_dir)) {
  dir.create(reports_dir, recursive = TRUE)
  cat("Created directory:", reports_dir, "\n")
}

# ---------- Table 1: Kruskal-Wallis results for WWTP_size ----------
table1 <- data.frame(
  Test = "CEC_loads vs. WWTP_size",
  Chi_squared = round(chi_sq1, 4),
  df = df1,
  p_value = round(p_val1, 6),
  Significance_at_alpha_0.05 = ifelse(significant1, "Yes", "No"),
  stringsAsFactors = FALSE
)

# Add post-hoc results if available
if (significant1 && !is.null(posthoc1)) {
  ph_table1 <- posthoc1$res
  ph_table1$Significant_adj <- ph_table1$P.adj < 0.05
  # Rename columns for clarity
  colnames(ph_table1)[colnames(ph_table1) == "Comparison"] <- "Pairwise_Comparison"
  colnames(ph_table1)[colnames(ph_table1) == "Z"] <- "Z_statistic"
  colnames(ph_table1)[colnames(ph_table1) == "P.unadj"] <- "Unadjusted_p"
  colnames(ph_table1)[colnames(ph_table1) == "P.adj"] <- "Bonferroni_adjusted_p"
  colnames(ph_table1)[colnames(ph_table1) == "Significant_adj"] <- "Significant_at_alpha_0.05"
}

# Save Table 1
table1_path <- file.path(tables_dir, "Table1_KW_WWTPsize.csv")
write.csv(table1, table1_path, row.names = FALSE)
cat("\nTable 1 saved to:", table1_path, "\n")

if (significant1 && !is.null(posthoc1)) {
  ph_table1_path <- file.path(tables_dir, "Table1_KW_WWTPsize_posthoc.csv")
  write.csv(ph_table1, ph_table1_path, row.names = FALSE)
  cat("Post-hoc table saved to:", ph_table1_path, "\n")
}

# ---------- Table 2: Kruskal-Wallis results for Community_size ----------
table2 <- data.frame(
  Test = "CEC_loads vs. Community_size",
  Chi_squared = round(chi_sq2, 4),
  df = df2,
  p_value = round(p_val2, 6),
  Significance_at_alpha_0.05 = ifelse(significant2, "Yes", "No"),
  stringsAsFactors = FALSE
)

# Add post-hoc results if available
if (significant2 && !is.null(posthoc2)) {
  ph_table2 <- posthoc2$res
  ph_table2$Significant_adj <- ph_table2$P.adj < 0.05
  colnames(ph_table2)[colnames(ph_table2) == "Comparison"] <- "Pairwise_Comparison"
  colnames(ph_table2)[colnames(ph_table2) == "Z"] <- "Z_statistic"
  colnames(ph_table2)[colnames(ph_table2) == "P.unadj"] <- "Unadjusted_p"
  colnames(ph_table2)[colnames(ph_table2) == "P.adj"] <- "Bonferroni_adjusted_p"
  colnames(ph_table2)[colnames(ph_table2) == "Significant_adj"] <- "Significant_at_alpha_0.05"
}

# Save Table 2
table2_path <- file.path(tables_dir, "Table2_KW_Communitysize.csv")
write.csv(table2, table2_path, row.names = FALSE)
cat("Table 2 saved to:", table2_path, "\n")

if (significant2 && !is.null(posthoc2)) {
  ph_table2_path <- file.path(tables_dir, "Table2_KW_Communitysize_posthoc.csv")
  write.csv(ph_table2, ph_table2_path, row.names = FALSE)
  cat("Post-hoc table saved to:", ph_table2_path, "\n")
}

# Save combined summary table
combined_path <- file.path(tables_dir, "Combined_Summary_WWTP_Community.csv")
write.csv(combined_summary, combined_path, row.names = FALSE)
cat("Combined summary table saved to:", combined_path, "\n")

# ==============================================================================
# 8. GENERATE PLAIN-TEXT REPORT FILES (.TXT FORMAT)
# ==============================================================================

# ---------- Report 1: WWTP_size ----------
report1_lines <- c(
  "=================================================================================",
  "          KRUSKAL-WALLIS TEST REPORT: CEC LOADS vs. WWTP SIZE",
  "=================================================================================",
  paste("Date and time of analysis:", Sys.time()),
  paste("Data file:", file_path),
  paste("Total observations (after cleaning):", nrow(df_clean)),
  "",
  "---------------------------------------------------------------------------------",
  "SUMMARY STATISTICS BY WWTP SIZE",
  "---------------------------------------------------------------------------------",
  paste(capture.output(print(summary_wwtp)), collapse = "\n"),
  "",
  "---------------------------------------------------------------------------------",
  "KRUSKAL-WALLIS TEST RESULTS",
  "---------------------------------------------------------------------------------",
  paste("Test name: Kruskal-Wallis rank sum test"),
  paste("Response variable: CEC_loads"),
  paste("Predictor variable: WWTP_size (levels:", 
        paste(levels(df_clean$WWTP_size), collapse = ", "), ")"),
  paste("Chi-squared statistic:", round(chi_sq1, 4)),
  paste("Degrees of freedom:", df1),
  paste("p-value:", format(p_val1, scientific = FALSE, digits = 6)),
  "",
  "---------------------------------------------------------------------------------",
  "INTERPRETATION",
  "---------------------------------------------------------------------------------"
)

if (significant1) {
  report1_lines <- c(report1_lines,
                     paste("There is a STATISTICALLY SIGNIFICANT difference in CEC loads across",
                           "WWTP size categories (χ² =", round(chi_sq1, 4), 
                           ", df =", df1, ", p =", format(p_val1, scientific = FALSE, digits = 6), ")."),
                     "",
                     "At α = 0.05, the null hypothesis of equal medians is REJECTED.",
                     "",
                     "---------------------------------------------------------------------------------",
                     "POST-HOC ANALYSIS (Bonferroni-Dunn pairwise comparisons)",
                     "---------------------------------------------------------------------------------"
  )
  
  if (!is.null(posthoc1)) {
    ph_text1 <- capture.output(print(ph_results1))
    report1_lines <- c(report1_lines, ph_text1, "")
    
    # Identify significant pairs
    sig_pairs1 <- ph_results1[ph_results1$P.adj < 0.05, ]
    if (nrow(sig_pairs1) > 0) {
      report1_lines <- c(report1_lines,
                         "SIGNIFICANT PAIRWISE DIFFERENCES (α = 0.05):",
                         paste(sig_pairs1$Comparison, " (adjusted p =", 
                               format(sig_pairs1$P.adj, scientific = FALSE, digits = 4), ")"),
                         ""
      )
    } else {
      report1_lines <- c(report1_lines,
                         "No pairwise comparisons were significant after Bonferroni correction.",
                         ""
      )
    }
    
    # Biological/environmental meaning
    report1_lines <- c(report1_lines,
                       "---------------------------------------------------------------------------------",
                       "ENVIRONMENTAL/BIOLOGICAL INTERPRETATION",
                       "---------------------------------------------------------------------------------",
                       "The significant overall difference suggests that WWTP size is associated with",
                       "different CEC loads in biosolids. This may reflect differences in treatment",
                       "technologies, sludge retention times, or influent characteristics between",
                       "small and large facilities. The specific pairwise differences (if any)",
                       "indicate which size classes are driving this effect.",
                       "",
                       "If larger plants show higher median loads, this could indicate:",
                       "  1) More complex influent mixtures in larger catchments",
                       "  2) Different sludge management practices",
                       "  3) Accumulation of persistent CECs over longer treatment chains",
                       "",
                       "Conversely, if smaller plants show higher loads, this might suggest:",
                       "  1) Less efficient treatment in small facilities",
                       "  2) Higher variability or shock loads in small systems",
                       "  3) Different sludge age or processing conditions"
    )
  }
  
} else {
  report1_lines <- c(report1_lines,
                     paste("There is NO STATISTICALLY SIGNIFICANT difference in CEC loads across",
                           "WWTP size categories (χ² =", round(chi_sq1, 4), 
                           ", df =", df1, ", p =", format(p_val1, scientific = FALSE, digits = 6), ")."),
                     "",
                     "At α = 0.05, we FAIL TO REJECT the null hypothesis of equal medians.",
                     "",
                     "---------------------------------------------------------------------------------",
                     "ENVIRONMENTAL/BIOLOGICAL INTERPRETATION",
                     "---------------------------------------------------------------------------------",
                     "The lack of significant difference suggests that WWTP size does not strongly",
                     "influence CEC loads in biosolids. This may indicate that treatment processes",
                     "are similar across size classes, or that variation within each size class",
                     "is larger than variation between classes.",
                     "",
                     "Other factors such as treatment type, sludge age, or seasonal variability",
                     "may be more important determinants of CEC loads than plant capacity alone."
  )
}

report1_lines <- c(report1_lines,
                   "",
                   "---------------------------------------------------------------------------------",
                   "END OF REPORT",
                   "================================================================================="
)

# Write Report 1
report1_path <- file.path(reports_dir, "Report1_KW_WWTPsize.txt")
writeLines(report1_lines, report1_path)
cat("Report 1 saved to:", report1_path, "\n")

# ---------- Report 2: Community_size ----------
report2_lines <- c(
  "=================================================================================",
  "          KRUSKAL-WALLIS TEST REPORT: CEC LOADS vs. COMMUNITY SIZE",
  "=================================================================================",
  paste("Date and time of analysis:", Sys.time()),
  paste("Data file:", file_path),
  paste("Total observations (after cleaning):", nrow(df_clean)),
  "",
  "---------------------------------------------------------------------------------",
  "SUMMARY STATISTICS BY COMMUNITY SIZE",
  "---------------------------------------------------------------------------------",
  paste(capture.output(print(summary_community)), collapse = "\n"),
  "",
  "---------------------------------------------------------------------------------",
  "KRUSKAL-WALLIS TEST RESULTS",
  "---------------------------------------------------------------------------------",
  paste("Test name: Kruskal-Wallis rank sum test"),
  paste("Response variable: CEC_loads"),
  paste("Predictor variable: Community_size (levels:", 
        paste(levels(df_clean$Community_size), collapse = ", "), ")"),
  paste("Chi-squared statistic:", round(chi_sq2, 4)),
  paste("Degrees of freedom:", df2),
  paste("p-value:", format(p_val2, scientific = FALSE, digits = 6)),
  "",
  "---------------------------------------------------------------------------------",
  "INTERPRETATION",
  "---------------------------------------------------------------------------------"
)

if (significant2) {
  report2_lines <- c(report2_lines,
                     paste("There is a STATISTICALLY SIGNIFICANT difference in CEC loads across",
                           "Community size categories (χ² =", round(chi_sq2, 4), 
                           ", df =", df2, ", p =", format(p_val2, scientific = FALSE, digits = 6), ")."),
                     "",
                     "At α = 0.05, the null hypothesis of equal medians is REJECTED.",
                     "",
                     "---------------------------------------------------------------------------------",
                     "POST-HOC ANALYSIS (Bonferroni-Dunn pairwise comparisons)",
                     "---------------------------------------------------------------------------------"
  )
  
  if (!is.null(posthoc2)) {
    ph_text2 <- capture.output(print(ph_results2))
    report2_lines <- c(report2_lines, ph_text2, "")
    
    # Identify significant pairs
    sig_pairs2 <- ph_results2[ph_results2$P.adj < 0.05, ]
    if (nrow(sig_pairs2) > 0) {
      report2_lines <- c(report2_lines,
                         "SIGNIFICANT PAIRWISE DIFFERENCES (α = 0.05):",
                         paste(sig_pairs2$Comparison, " (adjusted p =", 
                               format(sig_pairs2$P.adj, scientific = FALSE, digits = 4), ")"),
                         ""
      )
    } else {
      report2_lines <- c(report2_lines,
                         "No pairwise comparisons were significant after Bonferroni correction.",
                         ""
      )
    }
    
    # Biological/environmental meaning
    report2_lines <- c(report2_lines,
                       "---------------------------------------------------------------------------------",
                       "ENVIRONMENTAL/BIOLOGICAL INTERPRETATION",
                       "---------------------------------------------------------------------------------",
                       "The significant overall difference suggests that community size is associated",
                       "with different CEC loads in biosolids. This likely reflects variation in",
                       "contaminant sources, consumption patterns, and wastewater composition",
                       "between different population sizes.",
                       "",
                       "If larger communities show higher median loads, this could indicate:",
                       "  1) Higher per-capita pharmaceutical and personal care product use",
                       "  2) Greater diversity of CEC sources (hospitals, industry, agriculture)",
                       "  3) More complex sewer networks with longer transport times",
                       "",
                       "If smaller communities show higher loads, this might suggest:",
                       "  1) Less dilution from industrial or commercial sources",
                       "  2) Higher proportion of septic or decentralized systems",
                       "  3) Different demographic profiles (e.g., higher elderly population)"
    )
  }
  
} else {
  report2_lines <- c(report2_lines,
                     paste("There is NO STATISTICALLY SIGNIFICANT difference in CEC loads across",
                           "Community size categories (χ² =", round(chi_sq2, 4), 
                           ", df =", df2, ", p =", format(p_val2, scientific = FALSE, digits = 6), ")."),
                     "",
                     "At α = 0.05, we FAIL TO REJECT the null hypothesis of equal medians.",
                     "",
                     "---------------------------------------------------------------------------------",
                     "ENVIRONMENTAL/BIOLOGICAL INTERPRETATION",
                     "---------------------------------------------------------------------------------",
                     "The lack of significant difference suggests that community size does not",
                     "strongly influence CEC loads in biosolids. This may indicate that CEC",
                     "sources are relatively uniform across community sizes, or that treatment",
                     "processes effectively homogenize influent differences.",
                     "",
                     "Individual consumption patterns, industrial inputs, or seasonal factors",
                     "may be more important determinants than raw population size alone."
  )
}

report2_lines <- c(report2_lines,
                   "",
                   "---------------------------------------------------------------------------------",
                   "END OF REPORT",
                   "================================================================================="
)

# Write Report 2
report2_path <- file.path(reports_dir, "Report2_KW_Communitysize.txt")
writeLines(report2_lines, report2_path)
cat("Report 2 saved to:", report2_path, "\n")

# ==============================================================================
# 9. CONSOLE OUTPUT SUMMARY
# ==============================================================================

cat("\n")
cat("=================================================================================\n")
cat("                          ANALYSIS COMPLETE\n")
cat("=================================================================================\n")
cat("\nOUTPUT FILES GENERATED:\n")
cat("1. Tables (CSV):\n")
cat("   -", table1_path, "\n")
cat("   -", table2_path, "\n")
if (significant1) cat("   -", file.path(tables_dir, "Table1_KW_WWTPsize_posthoc.csv"), "\n")
if (significant2) cat("   -", file.path(tables_dir, "Table2_KW_Communitysize_posthoc.csv"), "\n")
cat("   -", combined_path, "\n")
cat("\n2. Reports (TXT):\n")
cat("   -", report1_path, "\n")
cat("   -", report2_path, "\n")
cat("\n=================================================================================\n")
cat("                    SCRIPT EXECUTED SUCCESSFULLY\n")
cat("=================================================================================\n")

# End of script
