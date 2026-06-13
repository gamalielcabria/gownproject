library(tidyverse)

# Create mock input dataset: 24 samples x 10 OTUs
set.seed(123)

otu_df <- tibble(
  SampleID = paste0("Sample_", sprintf("%02d", 1:24)),
  OTU_001 = rpois(24, lambda = 105),
  OTU_002 = rpois(24, lambda = 50),
  OTU_003 = rpois(24, lambda = 2),
  OTU_004 = rpois(24, lambda = 295),
  OTU_005 = rpois(24, lambda = 14),
  OTU_006 = rpois(24, lambda = 78),
  OTU_007 = rpois(24, lambda = 6),
  OTU_008 = rpois(24, lambda = 148),
  OTU_009 = rpois(24, lambda = 27),
  OTU_010 = rpois(24, lambda = 1)
)

otu_df

# Metadata:
# 2 seasons
# 2 locations per season
# 3 weeks per location-season
# 2 duplicate samples per week
metadata_df <- expand_grid(
  Season = c("Spring", "Fall"),
  Location = c("Site_A", "Site_B"),
  Week = paste0("Week_", 1:3),
  Duplicate = c("A", "B")
) %>%
  mutate(
    SampleID = paste0("Sample_", sprintf("%02d", row_number()))
  ) %>%
  select(SampleID, Season, Location, Week, Duplicate)

metadata_df

taxonomy_df <- tribble(
  ~OTU,      ~Phylum,
  "OTU_001", "Bacillota",
  "OTU_002", "Pseudomonadota",
  "OTU_003", "Bacteroidota",
  "OTU_004", "Actinomycetota",
  "OTU_005", "Bacillota",
  "OTU_006", "Chloroflexota",
  "OTU_007", "Pseudomonadota",
  "OTU_008", "Bacteroidota",
  "OTU_009", "Actinomycetota",
  "OTU_010", "Acidobacteriota"
)

taxonomy_df

# Combine into long DF
plot_df <- otu_df %>%
  pivot_longer(
    cols = starts_with("OTU_"),
    names_to = "OTU",
    values_to = "Count"
  ) %>%
  left_join(metadata_df, by = "SampleID") %>%
  left_join(taxonomy_df, by = "OTU") %>%
  group_by(SampleID) %>%
  mutate(RelAbund = Count / sum(Count)) %>%
  ungroup()

plot_df

# ========================================
# This is the old plot using ggh4x nested
# ========================================
# Plot in a nested abundance:
# library(tidyverse)
# library(ggh4x)
# library(grid)
# 
# p <- ggplot(
#   plot_df,
#   aes(x = Duplicate, y = RelAbund, fill = Phylum)
# ) +
#   geom_col(width = 0.95) +
#   
#   facet_nested(
#     . ~ Season + Location + Week,
#     space = "free_x",
#     scales = "free_x",
#     remove_labels = FALSE,
#     nest_line = element_line(
#       colour = c("grey50"),
#       linetype = c("dotted"),
#       linewidth = c(0.5)
#     ),
#     solo_line = TRUE,
#     strip = strip_nested(
#       clip = "off",
#       background_x = elem_list_rect(
#         fill = c("white", "white", "white"),
#         color = c("white", "white", "white"),
#         linewidth = c(1, 0.5, 0.5)
#       ),
#       text_x = elem_list_text(
#         color = c("black", "black", "black"),
#         face = c("bold", "plain", "plain"),
#         size = c(11, 10, 6)
#       ),
#       by_layer_x = TRUE
#     )
#   ) +
#   
#   theme_bw() +
#   theme(
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     #panel.spacing.x = unit(c(0.5, 3, 0.5, 0.5, 3, 0.5), "mm"),
#     panel.spacing.x = unit(0.5, "lines"),
#     
#     panel.border = element_rect(
#       color = "white",
#       fill = NA,
#       linewidth = 0.5
#     ),
#     
#     axis.title.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     
#     legend.position = "right"
#   ) +
#   
#   labs(
#     y = "Relative abundance",
#     fill = "Phylum"
#   )
# 
# p

# ========================================================
# This is using ggh4x and cowplot to combine season plots
# ========================================================

library(tidyverse)
library(ggh4x)
library(grid)
library(patchwork)
library(cowplot)

# =========================
# 1. Dynamic season order
# =========================
season_levels <- plot_df %>%
  distinct(Season) %>%
  pull(Season)

# =========================
# 2. Season plot function
# =========================
make_season_plot <- function(season_name) {
  
  plot_df_season <- plot_df %>%
    filter(Season == season_name)
  
  ggplot(
    plot_df_season,
    aes(x = Duplicate, y = RelAbund, fill = Phylum)
  ) +
    geom_col(width = 0.95) +
    
    facet_nested(
      . ~ Location + Week,
      space = "free_x",
      scales = "free_x",
      remove_labels = FALSE,
      nest_line = element_line(
        colour = "grey50",
        linetype = "dotted",
        linewidth = 0.5
      ),
      solo_line = TRUE,
      strip = strip_nested(
        clip = "off",
        background_x = elem_list_rect(
          fill = c("white", "white"),
          color = c("white", "white"),
          linewidth = c(0.8, 0.5)
        ),
        text_x = elem_list_text(
          color = c("black", "black"),
          face = c("bold", "plain"),
          size = c(10, 6)
        ),
        by_layer_x = TRUE
      )
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = c(0, 0)
    ) +
    
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.spacing.x = unit(0.5, "lines"),
      
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      
      axis.ticks.x = element_blank(),
      
      legend.position = "none",
      
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 12
      ),
      
      plot.background = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1
      ),
      
      plot.margin = margin(0, 5, 0, 5)
    ) +
    
    labs(
      title = season_name,
      fill = "Phylum"
    )
}

# =========================
# 3. Make season plots with adjustable spacer
# =========================
season_plots <- map(season_levels, make_season_plot)

season_gap <- 0.02   # smaller = less space; try 0.01, 0.02, 0.05

season_plots_spaced <- list()

for (i in seq_along(season_plots)) {
  season_plots_spaced[[length(season_plots_spaced) + 1]] <- season_plots[[i]]
  
  if (i < length(season_plots)) {
    season_plots_spaced[[length(season_plots_spaced) + 1]] <- plot_spacer()
  }
}

season_widths <- c(
  rep(c(1, season_gap), length(season_plots) - 1),
  1
)

season_panel <- wrap_plots(
  season_plots_spaced,
  ncol = length(season_plots_spaced),
  widths = season_widths
)

# =========================
# 4. Separate y-axis plot
# =========================
axis_plot <- ggplot(
  data.frame(x = 1, y = c(0, 1)),
  aes(x = x, y = y)
) +
  geom_blank() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.title.y = element_text(size = 11),
    axis.text.y = element_text(size = 9),
    axis.ticks.y = element_line(),
    
    plot.background = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  labs(y = "Relative abundance")

# =========================
# 5. Separate legend plot
# =========================
legend_source <- ggplot(
  plot_df,
  aes(x = Duplicate, y = RelAbund, fill = Phylum)
) +
  geom_col() +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  labs(fill = "Phylum")

legend_only <- cowplot::get_legend(legend_source)

# =========================
# 6. Combine Final plot
# =========================
p <- axis_plot + season_panel + wrap_elements(legend_only) +
  plot_layout(
    widths = c(0.01, 1, 0.18)
  )

p
