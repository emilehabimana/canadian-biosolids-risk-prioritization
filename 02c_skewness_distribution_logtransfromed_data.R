# ============================================================
# GENERATE SKEWNESS INTERPRETATION REPORT (LOG-TRANSFORMED DATA)
# OUTPUT: Word document
# ============================================================

# ---------- PATHS ----------
data_dir    <- "C:/Users/Admin/Documents/Data/Tables"
reports_dir <- "C:/Users/Admin/Documents/Data/Reports"

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

cat("✓ Log-transformed skewness report saved to:\n",
    
