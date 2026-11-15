# ============================================================
# PCA – EC CLASS DRIVERS OF BETWEEN‑MUNICIPALITY VARIATION
# FINAL VERSION (Opioids and Stimulants label fixed)
# ============================================================

# ---------- PATHS ----------
data_dir    <- "C:/Users/Admin/Documents/Article_2/Data/Data"
tables_dir  <- "C:/Users/Admin/Documents/Article_2/Data/Tables"
figures_dir <- "C:/Users/Admin/Documents/Article_2/Data/Figures"

dir.create(tables_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- LIBRARIES ----------
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(stringr)

# ---------- LOAD MUNICIPALITY‑LEVEL PCA DATA ----------
ec_pca <- read.csv(
  file.path(data_dir, "ec_robust_scaled_logtrans_MUN_PCA.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ec_mat <- as.matrix(ec_pca[, -1])
mode(ec_mat) <- "numeric"
rownames(ec_mat) <- ec_pca$Municipality_code

# ---------- RUN PCA ----------
pca <- prcomp(ec_mat, center = TRUE, scale. = FALSE)

scores <- as.data.frame(pca$x)
scores$Municipality <- rownames(scores)
scores$Municipality_short <- str_replace(scores$Municipality, "^MUN", "M")

# ---------- CLUSTER MUNICIPALITIES IN PCA SPACE ----------
hc <- hclust(dist(scores[, c("PC1", "PC2")]), method = "ward.D2")
scores$Cluster <- factor(cutree(hc, k = 4))

write.csv(
  scores,
  file.path(tables_dir, "PCA_Municipality_Scores.csv"),
  row.names = FALSE
)

# ---------- LOAD EC → CLASS TABLE ----------
ec_class <- read.csv(
  file.path(data_dir, "ec_robust_scaled_logtrans_MUN_PCA_code_class.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(ec_class)[2] <- "EC"
colnames(ec_class)[4] <- "Class_code"
colnames(ec_class)[5] <- "Class_display"

# ---------- EXPLICIT CLASS NAME RECODING (AUTHORITATIVE) ----------
class_map <- c(
  "Abx_Amx"     = "Antibiotics and antimicrobials",
  "Oth_Pharm"   = "Other pharmaceuticals",
  "Op_Stim"     = "Opioids and Stimulants",
  "Trz_Herb"    = "Herbicides",
  "Horm"        = "Hormones",
  "ASM_AED"     = "Anticonvulsants",
  "Neonics_Fip" = "Neonicotinoids and fipronil",
  "PCPs"        = "Personal care products",
  "Bisph"       = "Bisphenols",
  "Ind_Chem"    = "Industrial chemicals",
  "ASs"         = "Artificial sweeteners",
  "Oth_Pest"    = "Other pesticides"
)

# ---------- PCA LOADINGS ----------
loadings <- as.data.frame(pca$rotation)
loadings$EC <- rownames(loadings)

loadings_class <- loadings %>%
  left_join(ec_class[, c("EC", "Class_code")], by = "EC") %>%
  mutate(
    Class_full = recode(Class_code, !!!class_map),
    Class_full = ifelse(is.na(Class_full),
                        "Antibiotics and antimicrobials",
                        Class_full)
  )

# ---------- CLASS CONTRIBUTIONS ----------
class_contrib <- loadings_class %>%
  pivot_longer(cols = starts_with("PC"),
               names_to = "PC",
               values_to = "Loading") %>%
  filter(PC %in% c("PC1", "PC2")) %>%
  group_by(Class_full, PC) %>%
  summarise(
    Contribution = sum(abs(Loading)),
    .groups = "drop"
  )

write.csv(
  class_contrib,
  file.path(tables_dir, "PCA_EC_Class_Contributions.csv"),
  row.names = FALSE
)

# ---------- FORCE EXACT Y‑AXIS ORDER (TOP → BOTTOM) ----------
ec_class_order <- c(
  "Antibiotics and antimicrobials",
  "Other pharmaceuticals",
  "Opioids and Stimulants",
  "Herbicides",
  "Hormones",
  "Anticonvulsants",
  "Neonicotinoids and fipronil",
  "Personal care products",
  "Bisphenols",
  "Industrial chemicals",
  "Artificial sweeteners",
  "Other pesticides"
)

class_contrib$Class_full <- factor(
  class_contrib$Class_full,
  levels = rev(ec_class_order)
)

# ---------- FIGURE 1: PCA SCORES ----------
p_scores <- ggplot(
  scores,
  aes(PC1, PC2, color = Cluster, label = Municipality_short)
) +
  geom_point(size = 3) +
  geom_text(
    size = 3,
    vjust = -0.9,
    family = "Arial",
    color = "black"
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, color = "black"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  ) +
  labs(
    x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "%)"),
    color = "Cluster"
  )

ggsave(
  file.path(figures_dir, "PCA_Municipality_Scores.png"),
  p_scores,
  width = 120/25.4,
  height = 90/25.4,
  dpi = 300
)

# ---------- FIGURE 2: EC CLASS CONTRIBUTIONS ----------
p_class <- ggplot(
  class_contrib,
  aes(Class_full, Contribution, fill = PC)
) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_y_discrete(drop = FALSE) +
  scale_fill_manual(
    values = c("PC1" = "#fb0101", "PC2" = "#70ffee")
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, color = "black"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  ) +
  labs(
    x = "EC class",
    y = "Sum of absolute loadings",
    fill = "PCA"
  )

ggsave(
  file.path(figures_dir, "PCA_EC_Class_Contributions.png"),
  p_class,
  width = 90/25.4,
  height = 90/25.4,
  dpi = 300
)

# ---------- ENSURE BAR ORDER: PC1 FIRST, PC2 SECOND ----------
class_contrib$PC <- factor(class_contrib$PC, levels = c("PC1", "PC2"))

# ---------- EC CLASS CONTRIBUTION BAR CHART ----------
p_class <- ggplot(
  class_contrib,
  aes(Class_full, Contribution, fill = PC)
) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_y_discrete(drop = FALSE) +
  scale_fill_manual(
    values = c("PC2" = "#70ffee", "PC1" = "#fb0101")
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    axis.text  = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9, color = "black"),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8)
  ) +
  labs(
    x = "EC class",
    y = "Sum of absolute loadings",
    fill = "PCA"
  )

ggsave(
  file.path(figures_dir, "PCA_EC_Class_Contributions.png"),
  p_class,
  width = 90/25.4,
  height = 90/25.4,
  dpi = 300
)

cat("✓ PCA script completed with final class labeling and ordering\n")