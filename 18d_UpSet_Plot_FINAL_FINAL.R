############################################################
# 0. SETUP
############################################################

setwd("C:/Users/Admin/Documents/Article_2/Data/Data")

dir.create("C:/Users/Admin/Documents/Article_2/Data/Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("C:/Users/Admin/Documents/Article_2/Data/Reports", showWarnings = FALSE, recursive = TRUE)

library(tidyverse)
library(cowplot)
library(UpSetR)

############################################################
# 1. READ DATA
############################################################

merged <- read_csv("C:/Users/Admin/Documents/Article_2/Data/Data/merged_rankings_FINAL.csv")

merged <- merged %>%
  arrange(Cumulative_ranking) %>%
  mutate(
    Compound_names = factor(Compound_names, levels = Compound_names)
  )

############################################################
# 2. FUNCTION TO BUILD ONE CUSTOM UPSET‑STYLE FIGURE
############################################################

make_upset_plot <- function(df_subset, fig_name) {
  
  df_subset <- df_subset %>%
    arrange(Cumulative_ranking) %>%
    mutate(Compound_names = factor(Compound_names, levels = Compound_names))
  
  # Long format for dot matrix
  dot_df <- df_subset %>%
    select(Compound_names, Ranking_RQ_EU, Ranking_RQ_CA, Ranking_HS, Ranking_TU) %>%
    pivot_longer(
      cols = starts_with("Ranking_"),
      names_to = "Criterion",
      values_to = "Rank"
    ) %>%
    mutate(
      Criterion = factor(
        Criterion,
        levels = c(
          "Ranking_RQ_EU",
          "Ranking_RQ_CA",
          "Ranking_HS",
          "Ranking_TU"
        ),
        labels = c(
          "RQ (European framework)",
          "RQ (Canadian framework)",
          "Hazard score",
          "Toxic unit"
        )
      ),
      Color = case_when(
        Criterion == "RQ (European framework)" ~ "blue",
        Criterion == "RQ (Canadian framework)" ~ "#FF2400",
        Criterion == "Hazard score"            ~ "yellow",
        Criterion == "Toxic unit"              ~ "cyan2"
      )
    )
  
  # TOP BAR CHART — strict alignment
  p_top <- ggplot(df_subset, aes(x = Compound_names, y = Cumulative_ranking)) +
    geom_col(fill = "grey30", width = 0.9) +
    scale_x_discrete(expand = c(0,0)) +     # alignment fix
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x  = element_blank(),
      axis.title.x = element_blank(),
      plot.margin  = margin(5, 5, 0, 5)
    ) +
    ylab("Cumulative ranking")
  
  # BOTTOM DOT MATRIX — strict alignment
  p_bottom <- ggplot(dot_df, aes(x = Compound_names, y = Criterion)) +
    geom_point(aes(color = Color), size = 3) +
    scale_color_identity() +
    scale_x_discrete(expand = c(0,0)) +     # alignment fix
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x  = element_text(angle = 40, hjust = 1, vjust = 1,
                                  family = "Arial", colour = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin  = margin(0, 5, 5, 5)
    )
  
  final_plot <- plot_grid(
    p_top, p_bottom,
    ncol = 1,
    rel_heights = c(1.2, 1)
  )
  
  ggsave(
    filename = file.path("C:/Users/Admin/Documents/Article_2/Data/Figures",
                         paste0(fig_name, ".png")),
    plot     = final_plot,
    width    = 190, units = "mm", dpi = 300
  )
}

############################################################
# 3. SPLIT INTO FOUR SUBSETS AND GENERATE FIGURES
############################################################

subset1 <- merged %>% slice(1:30)
subset2 <- merged %>% slice(31:60)
subset3 <- merged %>% slice(61:85)
subset4 <- merged %>% slice(86:110)

make_upset_plot(subset1, "Figure_1_1to30")
make_upset_plot(subset2, "Figure_2_31to60")
make_upset_plot(subset3, "Figure_3_61to85")
make_upset_plot(subset4, "Figure_4_86to110")

############################################################
# 4. GENERATE TEXT REPORT
############################################################

report_text <- "
UpSet-style ranking analysis of 110 contaminants

This report summarizes the comparative ranking of 110 contaminants across four criteria:
- RQ (European framework)
- RQ (Canadian framework)
- Hazard score
- Toxic unit

The input table 'merged_rankings_FINAL' contains, for each compound, the four individual rankings
and a cumulative ranking defined as:

Cumulative_ranking = Ranking_RQ_EU + Ranking_RQ_CA + Ranking_HS + Ranking_TU

Compounds were ordered from best (lowest cumulative ranking) to worst.

UpSet-style visualization
For each compound, an UpSet-style layout was constructed with:
- A top bar chart showing the cumulative ranking.
- A bottom 4-row dot matrix indicating the ranking position in each criterion:
  Row 1: RQ (European framework)
  Row 2: RQ (Canadian framework)
  Row 3: Hazard score
  Row 4: Toxic unit

Compound names are displayed along the x-axis, rotated by 40 degrees in Arial font.

Four figures were generated:
- Figure 1: compounds ranked 1–30
- Figure 2: compounds ranked 31–60
- Figure 3: compounds ranked 61–85
- Figure 4: compounds ranked 86–110

These visualizations help identify:
- Contaminants jointly prioritized by European and Canadian RQ frameworks.
- Cases where hazard score or toxic unit strongly influence overall priority.
- Divergences between frameworks that may require further investigation.

Overall, the UpSet-style approach provides an integrated, visually intuitive summary of multi-criteria
rankings, supporting transparent prioritization of contaminants for monitoring and management.
"

writeLines(
  report_text,
  "C:/Users/Admin/Documents/Article_2/Data/Reports/UpSet_ranking_report_from_merged_FINAL.txt"
)
