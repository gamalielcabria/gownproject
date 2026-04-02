#!/usr/bin/env Rscript

# =========================================================
# Sample overlap across years
# Input: 00_sampleoverlaps.csv
# Output:
#   01_year_overlap_upset.png
#   02_year_pairwise_overlap_heatmap.png
#   03_year_pairwise_overlap_counts.csv
#   04_well_presence_absence.csv
#   05_well_presence_dotplot.png
# =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexUpset)
  library(ggplot2)
  library(tibble)
})

input_folder <- "/home/glbcabria/Workbench/gownproject/P0-GOWN/300_TimeSeries/"

# ----------------------------
# 1) Read data
# ----------------------------
infile <- paste0(input_folder, "00_sampleoverlaps.csv")

df <- read.csv(
  infile,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)

# ----------------------------
# 2) Convert wide -> long
# ----------------------------
df_long <- df %>%
  pivot_longer(
    cols = everything(),
    names_to = "Year",
    values_to = "Well"
  ) %>%
  mutate(
    Well = str_trim(Well),
    Well = na_if(Well, ""),
    Well = na_if(Well, "NA")
  ) %>%
  filter(!is.na(Well))

# ----------------------------
# Optional cleanup
# ----------------------------
# If "GW381sic" is actually meant to be "GW381", uncomment below:
# df_long <- df_long %>%
#   mutate(Well = recode(Well, "GW381sic" = "GW381"))

# Remove duplicate Year-Well combinations if any
df_long <- df_long %>%
  distinct(Year, Well)

# Keep years in intended order
year_order <- colnames(df)
df_long$Year <- factor(df_long$Year, levels = year_order)

# ----------------------------
# 3) Build presence/absence matrix
# ----------------------------
presence_df <- df_long %>%
  mutate(Present = 1L) %>%
  pivot_wider(
    names_from = Year,
    values_from = Present,
    values_fill = 0
  ) %>%
  arrange(Well)

# Save presence/absence table
write.csv(
  presence_df,
  paste0(input_folder, "04_well_presence_absence.csv"),
  row.names = FALSE
)

# ----------------------------
# 4) UpSet plot using ComplexUpset
# ----------------------------
# ComplexUpset expects one row per item and one binary column per set.
# We keep the Well column for reference, but only intersect over year columns.
# Convert to logical cleanly (IMPORTANT)

# Convert to logical cleanly
# ----------------------------
# 4) UpSet plot using ComplexUpset
# ----------------------------

presence_plot_df <- presence_df %>%
  mutate(across(all_of(year_order), ~ . == 1))

upset_plot <- ComplexUpset::upset(
  presence_plot_df,
  intersect = year_order,
  sort_intersections_by = c("degree", "cardinality"),
  sort_sets = FALSE,
  width_ratio = 0.25,
  height_ratio = 0.9,
  min_size = 1,

  matrix = intersection_matrix(
    geom = geom_point(size = 2.5),
    segment = geom_segment(linewidth = 0.5),
    outline_color = list(
      active   = "black",
      inactive = scales::alpha("grey60", 0.05)
    )
  ),

  base_annotations = list(
    "Intersection size" = intersection_size(
      text = list(size = 3)
    ) +
      ylab("Number of shared wells")
  ),

  set_sizes = upset_set_size() +
    ylab("Wells per year")
) +
  ggtitle("Overlap of wells across sampling years") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12)
  )

ggsave(
  filename = paste0(input_folder, "01_year_overlap_upset.png"),
  plot = upset_plot,
  width = 16,
  height = 7,
  dpi = 300,
  bg = "white"
)

# ----------------------------
# 5) Pairwise overlap counts
# ----------------------------
mat <- presence_df %>%
  column_to_rownames("Well") %>%
  as.matrix()

# crossproduct gives year-by-year shared well counts
overlap_mat <- t(mat) %*% mat

overlap_df <- as.data.frame(overlap_mat) %>%
  rownames_to_column("Year1")

write.csv(
  overlap_df,
  paste0(input_folder, "03_year_pairwise_overlap_counts.csv"),
  row.names = FALSE
)

# ----------------------------
# 6) Pairwise overlap heatmap
# ----------------------------
heatmap_df <- as.data.frame(as.table(overlap_mat)) %>%
  rename(
    Year1 = Var1,
    Year2 = Var2,
    Shared_Wells = Freq
  )

heatmap_plot <- ggplot(heatmap_df, aes(x = Year1, y = Year2, fill = Shared_Wells)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Shared_Wells), size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  coord_fixed() +
  labs(
    title = "Pairwise overlap of wells between years",
    x = "Year",
    y = "Year",
    fill = "Shared wells"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  filename = paste0(input_folder, "02_year_pairwise_overlap_heatmap.png"),
  plot = heatmap_plot,
  width = 8,
  height = 7,
  dpi = 300
)

# ----------------------------
# 7) Presence dot plot by well across years
# ----------------------------
presence_long <- presence_df %>%
  pivot_longer(
    cols = -Well,
    names_to = "Year",
    values_to = "Present"
  ) %>%
  mutate(
    Year = factor(Year, levels = year_order)
  )

# Order wells by total number of years present, then alphabetically
well_order <- presence_df %>%
  mutate(n_years_present = rowSums(across(-Well))) %>%
  arrange(desc(n_years_present), Well) %>%
  pull(Well)

presence_long <- presence_long %>%
  mutate(
    Well = factor(Well, levels = rev(well_order))
  )

presence_plot <- ggplot(
  subset(presence_long, Present == 1),
  aes(x = Year, y = Well)
) +
  geom_point(
    size = 4.5,
    shape = 22,
    fill = "steelblue",
    color = "black",
    stroke = 0.5
  ) +
  scale_x_discrete(position = "top") +
  labs(
    title = "Presence of wells across sampling years",
    x = "Year",
    y = "Well"
  ) +
  theme_minimal(base_size = 11) +
  theme(
      panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, margin = margin(b = 2)),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )

n_wells <- n_distinct(presence_long$Well)
plot_height <- max(8, n_wells * 0.18)

ggsave(
  filename = paste0(input_folder, "05_well_presence_dotplot.png"),
  plot = presence_plot,
  width = 6,
  height = plot_height,
  dpi = 300,
  limitsize = FALSE
)

# ----------------------------
# 8) Console summary
# ----------------------------
cat("Done.\n")
cat("Files written:\n")
cat("  01_year_overlap_upset.png\n")
cat("  02_year_pairwise_overlap_heatmap.png\n")
cat("  03_year_pairwise_overlap_counts.csv\n")
cat("  04_well_presence_absence.csv\n")
cat("  05_well_presence_dotplot.png\n")