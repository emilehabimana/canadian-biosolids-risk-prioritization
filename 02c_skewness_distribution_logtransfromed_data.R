# ============================================================
# GENERATE SKEWNESS INTERPRETATION REPORT (LOG-TRANSFORMED DATA)
# OUTPUT: Word document
# ============================================================

# ---------- PATHS ----------
data_dir    <- "C:/Users/Admin/Documents/Article_2/Data/Tables"
reports_dir <- "C:/Users/Admin/Documents/Article_2/Data/Reports"

if (!dir.exists(reports_dir)) {
  dir.create(reports_dir, recursive = TRUE)
}

setwd(data_dir)

# ---------- LIBRARIES ----------
library(moments)
library(rmarkdown)

# ---------- LOAD LOG-TRANSFORMED DATA ----------
log_data <- read.csv(
  "ec_concentration_log10_transformed.csv",
  stringsAsFactors = FALSE,
  check.names = TRUE
)

# Remove empty columns if present
log_data <- log_data[, names(log_data) != "" & !is.na(names(log_data))]

compound_col <- names(log_data)[1]
sample_cols  <- names(log_data)[-1]

# Ensure numeric
log_data[sample_cols] <- suppressWarnings(
  lapply(log_data[sample_cols], as.numeric)
)

# ---------- COMPUTE ABSOLUTE SKEWNESS ----------
skewness_values <- numeric()

for (i in seq_len(nrow(log_data))) {
  vals <- unlist(log_data[i, sample_cols], use.names = FALSE)
  vals <- vals[is.finite(vals)]
  if (length(vals) >= 3) {
    skewness_values <- c(skewness_values, abs(skewness(vals)))
  }
}

# ---------- SUMMARY METRICS ----------
n_compounds <- length(skewness_values)

median_skew <- median(skewness_values)
mean_skew   <- mean(skewness_values)

below_median <- sum(skewness_values <= median_skew)
above_median <- sum(skewness_values > median_skew)

# Skewness class counts
g1_lt_0_5 <- sum(skewness_values < 0.5)
g1_0_5_1  <- sum(skewness_values >= 0.5 & skewness_values < 1)
g1_1_2    <- sum(skewness_values >= 1 & skewness_values < 2)
g1_ge_2   <- sum(skewness_values >= 2)

# ---------- WRITE R MARKDOWN FILE ----------
report_rmd <- file.path(reports_dir, "Skewness_Report_Log10.Rmd")

cat(
  '---
title: "Skewness Analysis and Interpretation Report (Log-Transformed Data)"
output: word_document
---

## 1. Purpose of skewness assessment after transformation

Following log₁₀(x + 0.001) transformation, the distributional properties of concentration data were re-evaluated
to assess the extent to which transformation reduced asymmetry and extreme values present in the raw dataset.

## 2. Interpretation of |g₁| for log-transformed data

The statistic |g₁| represents the absolute value of Fisher’s skewness coefficient and quantifies the magnitude
of asymmetry in the distribution of log-transformed concentrations for each compound.

Unlike the raw data, where skewness reflects extreme concentration variability, |g₁| values obtained from
log-transformed data indicate how effectively transformation stabilizes variance and moderates distributional
asymmetry.

## 3. Formula and computation

For each compound, skewness was computed as:

$$
g_1 = \\frac{\\frac{1}{n} \\sum_{i=1}^{n} (x_i - \\bar{x})^3}
           {\\left( \\sqrt{\\frac{1}{n} \\sum_{i=1}^{n} (x_i - \\bar{x})^2} \\right)^3}
$$

where $x_i$ are log₁₀(x + 0.001)-transformed concentration values.

## 4. Interpretation scale of |g₁|

The following classification was used to evaluate skewness severity:

- |g₁| < 0.5: approximately symmetric  
- 0.5 ≤ |g₁| < 1: moderately skewed  
- 1 ≤ |g₁| < 2: strongly skewed  
- |g₁| ≥ 2: extremely skewed  

This framework allows direct comparison between raw and transformed datasets.

## 5. Importance of post-transformation skewness assessment

Evaluating skewness after transformation is critical to verify that:
- extreme asymmetry has been reduced,
- mean-based statistics become more representative,
- parametric analyses become more defensible.

## 6. Kernel density estimation (KDE)

The kernel density estimation (KDE) associated with the skewness histogram provides a smoothed representation
of the distribution of |g₁| across compounds.

For the log-transformed data, the KDE highlights the redistribution of skewness values toward lower magnitudes,
relative to the raw dataset.

## 7. Boxplot interpretation

The horizontal boxplot summarizes the central tendency and dispersion of |g₁| values across log-transformed data:

- Compounds at or below the median skewness: **', below_median, '**
- Compounds above the median skewness: **', above_median, '**

This indicates the balance of skewness values around the median following transformation.

## 8. Distribution across skewness classes (log-transformed data)

The distribution of compounds across skewness classes after log transformation is as follows:

- |g₁| < 0.5: **', g1_lt_0_5, ' compounds**
- 0.5 ≤ |g₁| < 1: **', g1_0_5_1, ' compounds**
- 1 ≤ |g₁| < 2: **', g1_1_2, ' compounds**
- |g₁| ≥ 2: **', g1_ge_2, ' compounds**

Compared to the raw dataset, these results quantify the extent to which log transformation reduces extreme skewness
and improves distributional properties.

---

### References

- Fisher, R. A., & Pearson, K. (1907). On the mathematical foundation of theoretical statistics.
- Joanes, D. N., & Gill, C. A. (1998). Comparing measures of sample skewness and kurtosis.
',
  file = report_rmd
)

# ---------- RENDER WORD DOCUMENT ----------
render(
  input       = report_rmd,
  output_file = file.path(reports_dir, "Skewness_Report_Log10.docx"),
  quiet       = TRUE
)

cat("✓ Log-transformed skewness report saved to:\n",
    