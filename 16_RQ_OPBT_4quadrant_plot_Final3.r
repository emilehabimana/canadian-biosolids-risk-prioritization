############################################################
# OPBT × RQ Plot with Custom Piecewise Axis (FINAL CLEAN)
############################################################

# -------------------------
# Load required libraries
# -------------------------
library(ggplot2)
library(readr)
library(dplyr)
library(ggrepel)
library(scales)

# -------------------------
# Define paths
# -------------------------
data_path   <- "C:/Users/Admin/Documents/Data/Data/rq_opbt.csv"
figure_path <- "C:/Users/Admin/Documents/Data/Figures"

dir.create(figure_path, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Read data
# -------------------------
df <- read_csv(data_path)

df <- df %>%
  rename(
    compound_name = 1,
    compound_code = 2,
    RQ            = 3,
    OPBT          = 4
  )

# -------------------------
# Correct segment widths (mm → normalized)
# TOTAL = 190 mm
# -------------------------
w1 <- 60/190   # 0–0.001
w2 <- 35/190   # 0.001–0.01
w3 <- 45/190   # 0.01–0.1
w4 <- 25/190   # 0.1–1
w5 <- 25/190   # 1–45

# -------------------------
# Piecewise transformation
# -------------------------
transform_rq <- function(x) {
  case_when(
    x <= 0.001 ~ (x / 0.001) * w1,
    
    x <= 0.01  ~ w1 + (x - 0.001)/(0.009) * w2,
    
    x <= 0.1   ~ w1 + w2 + (x - 0.01)/(0.09) * w3,
    
    x <= 1     ~ w1 + w2 + w3 + (x - 0.1)/(0.9) * w4,
    
    TRUE       ~ w1 + w2 + w3 + w4 + (x - 1)/(45) * w5
  )
}

df$RQ_t <- transform_rq(df$RQ)

# -------------------------
# Segment boundaries
# -------------------------
b1 <- 0
b2 <- w1
b3 <- w1 + w2
b4 <- w1 + w2 + w3
b5 <- w1 + w2 + w3 + w4
b6 <- 1

# -------------------------
# Split dataset
# -------------------------
df_seg1  <- df %>% filter(RQ <= 0.001)
df_other <- df %>% filter(RQ > 0.001)

# -------------------------
# Point size scaling
# -------------------------
df$point_size <- rescale(df$RQ + df$OPBT, to = c(2.5, 6))

# -------------------------
# Modified labeling for segment 1 (0–0.001)
# Labels are nudged to the right and placed in empty horizontal space
# -------------------------
# Calculate available width within segment 1 for label placement
segment1_width <- b2 - b1
# Nudge factor: move labels to the right by 30-40% of segment width
nudge_factor <- segment1_width * 0.35

# For segment 1: create slightly nudged positions to the right
df_seg1 <- df_seg1 %>%
  mutate(
    # Move the label positions to the right within the same segment
    label_x = RQ_t + nudge_factor,
    # Ensure labels don't cross the 0.001 boundary (b2)
    label_x = pmin(label_x, b2 - 0.005)
  )

# -------------------------
# For all points: move name_codes next to dots (not on top)
# Create nudged positions for all data points
# -------------------------
# Calculate nudge based on RQ_t position and OPBT value
# For points with higher OPBT, nudge upward; for lower OPBT, nudge downward
# Also add horizontal nudge to the right for all points
df <- df %>%
  mutate(
    # Base horizontal nudge (right)
    nudge_x = 0.008,
    # Vertical nudge based on OPBT to avoid overlapping dots
    nudge_y = case_when(
      OPBT > 0.7 ~ -0.03,
      OPBT > 0.4 ~ 0.02,
      OPBT > 0.2 ~ -0.02,
      TRUE ~ 0.03
    ),
    # For segment 1, use the modified label positions
    label_x_final = ifelse(RQ <= 0.001, 
                           pmin(RQ_t + nudge_factor, b2 - 0.005),
                           RQ_t + nudge_x),
    label_y_final = OPBT + nudge_y
  )

# For other segments, prepare data for geom_text_repel with adjusted settings
df_other <- df %>% filter(RQ > 0.001) %>%
  mutate(
    label_x = RQ_t + nudge_x,
    label_y = OPBT + nudge_y
  )

# Update df_seg1 with final label positions
df_seg1 <- df %>% filter(RQ <= 0.001) %>%
  mutate(
    label_x = pmin(RQ_t + nudge_factor, b2 - 0.005),
    label_y = OPBT + nudge_y
  )

# -------------------------
# Build plot
# -------------------------
p <- ggplot() +
  
  # Horizontal dotted lines at y = 0.5 and y = 0.75 (full width)
  geom_hline(
    yintercept = 0.5,
    linetype = "dotted",
    color = "black",
    linewidth = 0.8,
    alpha = 0.7
  ) +
  
  geom_hline(
    yintercept = 0.75,
    linetype = "dotted",
    color = "black",
    linewidth = 0.8,
    alpha = 0.7
  ) +
  
  # Points
  geom_point(
    data = df,
    aes(x = RQ_t, y = OPBT, size = point_size),
    color = "#1A05A2",
    alpha = 0.85
  ) +
  
  # Labels for segments > 0.001 (with slight nudges to avoid dots)
  geom_text_repel(
    data = df_other,
    aes(x = RQ_t, y = OPBT, label = compound_code),
    size = 3,
    nudge_x = 0.008,
    nudge_y = 0.02,
    segment.color = "grey70",
    segment.size = 0.3,
    box.padding = 0.4,
    point.padding = 0.3,
    min.segment.length = 0.3,
    max.overlaps = 15
  ) +
  
  # Labels for segment 1 (WITH lines, extended to the right, within 0-0.001)
  geom_text_repel(
    data = df_seg1,
    aes(x = RQ_t, y = OPBT, label = compound_code),
    size = 3,
    # Force labels to be placed at the calculated right-shifted positions
    nudge_x = nudge_factor,
    direction = "both",
    segment.color = "grey40",
    segment.size = 0.3,
    min.segment.length = 0,
    xlim = c(b1 + segment1_width * 0.2, b2 - 0.005),
    ylim = c(0, 1),
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    force = 2
  ) +
  
  # -------------------------
  # COLORED AXIS SEGMENTS (USING EXACT COLOR NAMES)
  # Segment 1 (0–0.001): chartreuse3
  # Segment 2 (0.001–0.01): chartreuse2
  # Segment 3 (0.01–0.1): darkolivegreen1
  # Segment 4 (0.1–1): yellow
  # Segment 5 (1–40): firebrick1
  # -------------------------
  
  # Main colored axis segments (linewidth = 4 for better visibility)
  annotate("segment", x = b1, xend = b2, y = -0.04, yend = -0.04,
           color = "chartreuse3", linewidth = 4) +
  
  annotate("segment", x = b2, xend = b3, y = -0.04, yend = -0.04,
           color = "chartreuse2", linewidth = 4) +
  
  annotate("segment", x = b3, xend = b4, y = -0.04, yend = -0.04,
           color = "darkolivegreen1", linewidth = 4) +
  
  annotate("segment", x = b4, xend = b5, y = -0.04, yend = -0.04,
           color = "yellow", linewidth = 4) +
  
  annotate("segment", x = b5, xend = b6, y = -0.04, yend = -0.04,
           color = "firebrick1", linewidth = 4) +
  
  # Subtle background rectangles with alpha transparency to enhance color visibility
  annotate("rect", xmin = b1, xmax = b2, ymin = -0.045, ymax = -0.035,
           fill = "chartreuse3", alpha = 0.3) +
  annotate("rect", xmin = b2, xmax = b3, ymin = -0.045, ymax = -0.035,
           fill = "chartreuse2", alpha = 0.3) +
  annotate("rect", xmin = b3, xmax = b4, ymin = -0.045, ymax = -0.035,
           fill = "darkolivegreen1", alpha = 0.3) +
  annotate("rect", xmin = b4, xmax = b5, ymin = -0.045, ymax = -0.035,
           fill = "yellow", alpha = 0.3) +
  annotate("rect", xmin = b5, xmax = b6, ymin = -0.045, ymax = -0.035,
           fill = "firebrick1", alpha = 0.3) +
  
  # -------------------------
  # Segment separators
  # -------------------------
  geom_vline(
    xintercept = c(b2, b3, b4, b5),
    linetype = "dashed",
    color = "black",
    alpha = 0.6
  ) +
  
  # -------------------------
  # Custom axes
  # -------------------------
  scale_x_continuous(
    limits = c(0, 1),
    expand = c(0, 0),
    breaks = c(b1, b2, b3, b4, b5, b6),
    labels = c("0", "0.001", "0.01", "0.1", "1", "40")
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  
  scale_size_identity() +
  
  # -------------------------
  # Theme
  # -------------------------
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    # Add slight margin at bottom for colored axis segments
    plot.margin = margin(t = 5, r = 5, b = 10, l = 5, unit = "mm")
  ) +
  
  # No title
  labs(
    x = "Risk Quotient (RQ)",
    y = "OPBT Hazard Relevance Score"
  )

# -------------------------
# Save figure
# -------------------------
ggsave(
  filename = file.path(figure_path, "RQ_OPBT_Final.png"),
  plot     = p,
  width    = 190,
  height   = 190,
  units    = "mm",
  dpi      = 300
)

# Also save a PDF version for better quality if needed
ggsave(
  filename = file.path(figure_path, "RQ_OPBT_Final.pdf"),
  plot     = p,
  width    = 190,
  height   = 190,
  units    = "mm",
  dpi      = 300
)

############################################################
# END OF SCRIPT
############################################################
