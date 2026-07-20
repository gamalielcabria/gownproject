# test_nested_dotplot.R
# ---------------------------------------------------------------------------
# Quick, self-contained test harness for nested_dotplot().
#
# Expected files next to this script (or adjust the paths below):
#   R/nested_dotplot.R   - the function
#   mock_data.csv        - the wide-format mock data
#
# The mock data is WIDE (one column per ASV). nested_dotplot() wants LONG
# data (one row per Sample x ASV), so we pivot first.
#
# This version demonstrates:
#   * a different point colour for each ASV
#   * a different strip fill for each Type   (cDNA vs DNA)  -> strip_by_layer = FALSE
#   * a different title colour for each Category            -> title_colour = named vector
# ---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(patchwork)
library(cowplot)

# --- locate and source the function ----------------------------------------
fn_path <- if (file.exists("/home/glbcabria/Workbench/plotutilityscript/R/nested_dotplot.R")) {
  "/home/glbcabria/Workbench/plotutilityscript/R/nested_dotplot.R"
} else if (file.exists("nested_dotplot.R")) {
  "nested_dotplot.R"
} else {
  stop("Could not find nested_dotplot.R. Set the path manually.")
}
# If you have installed the package instead, replace the source() line with:
#   library(plotutilityscript)
source(fn_path)

# --- read the wide mock data ------------------------------------------------
csv_path <- if (file.exists("/home/glbcabria/Workbench/gownproject/Training/dotplot_test_data.csv")) {
  "/home/glbcabria/Workbench/gownproject/Training/dotplot_test_data.csv"
} else if (file.exists("data-raw/mock_data.csv")) {
  "data-raw/mock_data.csv"
} else {
  stop("Could not find mock_data.csv. Set the path manually.")
}
wide <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE) %>%
    mutate(Replicate = gsub("^", "R", Replicate))

# --- reshape to long --------------------------------------------------------
long <- wide |>
  tidyr::pivot_longer(
    cols = dplyr::starts_with("ASV"),
    names_to = "ASV",
    values_to = "RelAbund"
  ) |>
  dplyr::mutate(
    # x-axis and facet variables must be discrete
    Replicate = factor(Replicate, levels = c("R1", "R2", "R3")),
    Type      = factor(Type, levels = c("cDNA", "DNA")),
    Category  = factor(Category, levels = c("Diesel", "Gas99", "Gas95")),
    ASV       = factor(ASV, levels = paste0("ASV", 1:8))
  )

str(long)

# --- colour keys ------------------------------------------------------------
# One colour per ASV (used for the point colour aesthetic).
asv_cols <- c(
  ASV1 = "#1b9e77", ASV2 = "#d95f02", ASV3 = "#7570b3", ASV4 = "#e7298a",
  ASV5 = "#66a61e", ASV6 = "#e6ab02", ASV7 = "#a6761d", ASV8 = "#666666"
)

# --- the plot ---------------------------------------------------------------

p2 <- nested_dotplot(
  data              = long,
  x_col             = Replicate,
  y_col             = ASV,
  size_col          = RelAbund,
  
  # Orientation of faceted plots
  orientation       = "vertical",

  # Color parameters of inputs
  colour_col        = ASV,                     # <- different colour per ASV
  point_border      = "black",              # Color of the ring around the dots
  point_stroke      = 0.8,                      # Thickness of the ring around the dots
  split_col         = Type,                    # <- Type is now the title
  nested_cols       = c("Category"),           # <- Category is now the strip
  palette           = asv_cols,
  
  #Legend Labels
  colour_label      = "ASV",
  size_label        = "Relative abundance",
  legend_position   = "right",                  # <- legend under the plot
  
  #Axis text size
  x_text_size       = 12,
  y_text_size       = 20,

  # Spacing of y-axis, main panel, and legend
  layout_widths = c(0.01,1,0.18),

  # Rank of ASVs top to bottom
  order_by          = "abundance",

  # Strip Layers
  ## Strip Colors
  strip_by_layer    = FALSE,
  strip_colour      = c("blue","black","gray20"),             # Border Color for strips
  strip_fill        = c("#d8b365", "#5ab4ac", "#c2a5cf"),     # per-Category strip colours
  strip_text_colour = c("#3d2a05", "#053b37", "#2a0d33"),     # Text in strip color (can be just 1 color: "white")
  strip_text_face   = c("bold", "bold", "bold"),
  strip_text_size   = c(25,20,15,10),                               # Size of Strip Levels (e.g. Type Col, Category Col)
  
  ## Title is the First Level Name 
  title_colour      = c(cDNA = "#b35806", DNA = "#542788"),     # per-Type title colours
  title_position    = "right",                                      # Box Tiotle Position
  title_angle       = 270,                                          # <- ...rotated vertical

  # Themes                                                          # Takes in GGPLOT2 Themes
  plot_theme = ggplot2::theme(
        panel.grid.major.x = ggplot2::element_line(colour = "grey90", linewidth = 0.4)
    )
)

p2

ggsave2("/home/glbcabria/Workbench/plotutilityscript/man/figures/nested_dotplot-example-2.png", plot = p2, width = 12, height = 10, dpi = 300)
