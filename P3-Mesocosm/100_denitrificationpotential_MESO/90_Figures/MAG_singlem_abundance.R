# MAG_singlem_abundance.R
# This script matches the abundance of gtdbtk identified MAGs to the singleM relative abundance

# Load libraries
invisible(suppressWarnings(suppressPackageStartupMessages({
      library(tidyverse)
      library(scales)
      library(microeco)
      library(fuzzyjoin)
      library(janitor)
      library(ggbreak)
    })))

# Set input and output paths
assembly_input_path <- "gownproject/P3-Mesocosm/100_denitrificationpotential_MESO/20_assembly/04_gtdbtk/MLA_25MESO_SHMG_gtdbtk.bac120.summary-comebin.tsv"
singlem_input_path <- "gownproject/P3-Mesocosm/100_denitrificationpotential_MESO/30_taxonomicabundance/03_MESO25/03_MESO25_combined-profiles.summarise.wextras.tsv"
checkm_input_path <- "gownproject/P3-Mesocosm/100_denitrificationpotential_MESO/20_assembly/03_refinedbins/03_checkm/MLA_25MESO_SHMG_checkm_summary-comebin.tsv"
output_path <- "gownproject/P3-Mesocosm/100_denitrificationpotential_MESO/90_Figures/"

# FUNCTIONS
clean_rank <- function(x, prefix) {
  x <- str_trim(x)
  stripped <- str_remove(x, paste0("^", prefix))
  if_else(stripped == "", NA_character_, x)
}

# Load data
checkm_data <- read_tsv(checkm_input_path) %>%
  clean_names() %>%
  select(bin_id, completeness, contamination, strain_heterogeneity , genome_size_bp, number_contigs, n50_scaffolds, gc)

assembly_data <- read_tsv(assembly_input_path) %>%
  select(user_genome, classification, closest_genome_reference, closest_genome_ani) %>%
  separate(
    classification,
    into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = ";",
    fill = "right"
  ) %>%
  mutate(
    Domain  = clean_rank(Domain,  "d__"),
    Phylum  = clean_rank(Phylum,  "p__"),
    Class   = clean_rank(Class,   "c__"),
    Order   = clean_rank(Order,   "o__"),
    Family  = clean_rank(Family,  "f__"),
    Genus   = clean_rank(Genus,   "g__"),
    Species = clean_rank(Species, "s__")
  ) %>%
  unite(
    "Taxa",
    Domain, Phylum, Class, Order, Family, Genus, Species,
    sep = ";",
    na.rm = TRUE,
    remove = TRUE
  )

singlem_data <- read_tsv(singlem_input_path) %>%
    separate(taxonomy, into = c("Root","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "; ")%>%
    unite("Taxa", Domain, Phylum, Class, Order, Family, Genus, Species, sep = ";", na.rm = TRUE) %>%
    select(-Root)

# COMBINE data from assembly and singleM
assembly_singlem_data <- assembly_data %>%
  left_join(singlem_data, by = "Taxa") %>%
  mutate( 
    terminal_taxid = str_extract(Taxa, "[^;]+$"),
    parent_taxa = str_remove(Taxa, ";[^;]+$")
    )

# LVL0
assembly_singlem_data_wid <- assembly_singlem_data %>%
  filter(!is.na(sample))

# LVL1 
assembly_singlem_data_noid_lvl1 <- assembly_singlem_data %>%
  filter(is.na(sample)) %>%
  select(user_genome, closest_genome_reference, closest_genome_ani, terminal_taxid, parent_taxa) %>%
  mutate(Taxa = parent_taxa) %>%
  left_join(singlem_data, by = "Taxa")

# LVL2
left_tbl <- assembly_singlem_data_noid_lvl1 %>%
  filter(is.na(sample)) %>%
  select(user_genome, closest_genome_reference, closest_genome_ani, terminal_taxid, parent_taxa) %>%
  mutate(join_taxon = str_extract(parent_taxa, "[^;]+$"))

match_tbl <- left_tbl %>%
  mutate(row_id = row_number()) %>%
  transmute(
    row_id,
    join_taxon,
    match_data = map(
      join_taxon,
      ~ {
        if (is.na(.x) || .x == "") {
          singlem_data[0, ]
        } else {
          singlem_data %>%
            filter(str_detect(Taxa, fixed(.x)))
        }
      }
    )
  ) %>%
  unnest(match_data, keep_empty = TRUE)

assembly_singlem_data_noid_lvl2 <- left_tbl %>%
  mutate(row_id = row_number()) %>%
  left_join(match_tbl, by = c("row_id", "join_taxon"))

# LVL3
left_tbl <- assembly_singlem_data_noid_lvl2 %>%
  filter(is.na(sample)) %>%
  select(user_genome, closest_genome_reference, closest_genome_ani, terminal_taxid, parent_taxa) %>%
  mutate(join_taxon = terminal_taxid)

match_tbl <- left_tbl %>%
  mutate(row_id = row_number()) %>%
  transmute(
    row_id,
    join_taxon,
    match_data = map(
      join_taxon,
      ~ {
        if (is.na(.x) || .x == "") {
          singlem_data[0, ]
        } else {
          singlem_data %>%
            filter(str_detect(Taxa, fixed(.x)))
        }
      }
    )
  ) %>%
  unnest(match_data, keep_empty = TRUE)

assembly_singlem_data_noid_lvl3 <- left_tbl %>%
  mutate(row_id = row_number()) %>%
  left_join(match_tbl, by = c("row_id", "join_taxon"))

# combine all levels
assembly_singlem_data_combined <- bind_rows(
  assembly_singlem_data_wid %>% select(user_genome, closest_genome_reference, closest_genome_ani, terminal_taxid, parent_taxa,sample, relative_abundance, level),
  assembly_singlem_data_noid_lvl1 %>% select(user_genome, closest_genome_reference, closest_genome_ani, terminal_taxid, parent_taxa,sample, relative_abundance, level),
  assembly_singlem_data_noid_lvl2 %>% select(user_genome, closest_genome_reference, closest_genome_ani, terminal_taxid, parent_taxa,sample, relative_abundance, level),
  assembly_singlem_data_noid_lvl3 %>% select(user_genome, closest_genome_reference, closest_genome_ani, terminal_taxid, parent_taxa,sample, relative_abundance, level)
) %>%
filter(!is.na(sample)) %>%
  mutate(
  code = str_extract(sample, "(CO|SH|SS|SA)[0-9]+"),
  Substrate = case_when(
      str_starts(code, "CO") ~ "Coal",
      str_starts(code, "SH") ~ "Shale",
      str_starts(code, "SS") ~ "Sandstone",
      str_starts(code, "SA") ~ "Sand",
      TRUE ~ NA_character_
  )
) %>%
mutate(
  MAG_ID = paste0(
      "MAG_M",
      sprintf("%02d", as.integer(factor(user_genome)))
  )
) %>%
filter(
  !(terminal_taxid == 'g__JAEUWI01' & level == 'family')
) %>%
filter(
  !(terminal_taxid == 'f__R501' & level %in% c('family','genus'))
) %>%
arrange(MAG_ID) %>%
select(MAG_ID, code, Substrate, parent_taxa, terminal_taxid, relative_abundance, closest_genome_reference, closest_genome_ani, user_genome)

checkm_assembly_data <- checkm_data %>%
  left_join(
    assembly_singlem_data_combined %>% 
    group_by(MAG_ID, user_genome, parent_taxa, terminal_taxid) %>%
    summarise(max_relative_abundance = max(relative_abundance, na.rm = TRUE), .groups = "drop"),
    by = c("bin_id" = "user_genome")
  ) %>%
  unite(
    "Taxa",
    parent_taxa, terminal_taxid,
    sep = ";",
    na.rm = TRUE,
    remove = TRUE
  ) %>%
  select(MAG_ID, everything())

# ========================================
# SAVE TABLES
# ========================================
# write.csv(assembly_singlem_data_combined, file.path(output_path, "assembly_singlem_data_combined.tsv"))
# write.csv(checkm_assembly_data, file.path(output_path, "assembly_checkm_data.tsv"))

# ========================================
# PLOT TABLE
# ========================================
## Set Color Palette (Universal)
substrate_colors <- c(
  "Sandstone" = "#868686FF",      # grey (dark rock)
  "Coal" = "#0073C2FF",     # blue (cool sediment)
  "Shale" = "#CD534CFF", # red (iron-rich)
  "Sand" = "#EFC000FF"       # yellow (sand)
)

## Create a plot table by joining the combined assembly and singleM data with the checkm data
plot_table <- assembly_singlem_data_combined %>%
  select(-parent_taxa) %>%
  left_join(
    checkm_assembly_data %>% select(MAG_ID, completeness, contamination, genome_size_bp, number_contigs, n50_scaffolds, gc, Taxa),
    by = "MAG_ID"
  ) %>%
  mutate(Order = str_extract(Taxa, "o__[^;]+")) %>%
  rename(TaxID = terminal_taxid)

## PLOT IT!#####
plot_table <- plot_table %>%
  mutate(
    # x_label = paste0("(", gsub("o__", "", Order), ") ", MAG_ID)
    x_label = paste0("(", gsub("^[a-z]__", "", TaxID), ") ", MAG_ID)
  ) %>%
  arrange(Order) %>%  # ensure row order follows Order
  mutate(
    x_label = factor(x_label, levels = unique(x_label))
  )

# plot_mag_abundance <- ggplot(
#   plot_table,
#   aes(x = x_label, y = relative_abundance)
# ) +
#   geom_boxplot(
#     aes(fill = Substrate),
#     outlier.shape = NA,
#     alpha = 0.6
#   ) +
#   geom_jitter(
#     aes(color = Substrate),
#     width = 0.2,
#     alpha = 1
#   ) +
#   scale_fill_manual(values = substrate_colors) +
#   scale_color_manual(values = substrate_colors) +
#   facet_wrap(~Substrate, ncol = 1) +
#   labs(
#     x = "MAG ID (Order)",
#     y = "Relative Abundance (%)",
#     title = "Relative Abundance of MAGs across Substrates"
#   ) +
#   theme_bw() +
#   theme(
#     # aspect.ratio = 0.15,
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#     legend.position = "right"
#   )

# plot_mag_abundance

# ggsave(
#   plot = plot_mag_abundance,
#   filename = file.path(output_path, "MAG_singlem_abundance.png"),
#   width = 15,
#   height = 8,
#   dpi = 300
# )


# ===========================================
# PLOT WITH ORDER LABELS AND SEPARATION LINES
#===========================================
library(dplyr)
library(ggplot2)
library(patchwork)

#----------------------------
# 1) Prepare data
#----------------------------
plot_table2 <- plot_table %>%
  mutate(
    Order_clean = gsub("^[a-z]__", "", Order),
    TaxID_clean = gsub(" .*$", "", TaxID),
    x_label = paste0("(", TaxID_clean, ") ", MAG_ID)
    # if you want simpler labels, use:
    # x_label = MAG_ID
  ) %>%
  arrange(Order_clean, MAG_ID) %>%
  mutate(
    x_label = factor(x_label, levels = unique(x_label)),
    x_num   = as.numeric(x_label)
  )

#----------------------------
# 2) Data for order labels and boundaries
#----------------------------
order_centers <- plot_table2 %>%
  distinct(Order_clean, x_label, x_num) %>%
  group_by(Order_clean) %>%
  summarise(
    xmin   = min(x_num),
    xmax   = max(x_num),
    xmid   = mean(range(x_num)),
    .groups = "drop"
  )

order_boundaries <- order_centers %>%
  mutate(boundary = xmax + 0.5)

#----------------------------
# 3) Main plot
#----------------------------
p_main <- ggplot(
  plot_table2,
  aes(x = x_label, y = relative_abundance)
) +
  geom_boxplot(
    aes(fill = Substrate),
    outlier.shape = NA,
    alpha = 0.6
  ) +
  geom_jitter(
    aes(color = Substrate),
    width = 0.2,
    alpha = 1
  ) +
  # geom_vline(
  #   data = order_boundaries %>% filter(boundary < max(plot_table2$x_num) + 0.5),
  #   aes(xintercept = boundary),
  #   inherit.aes = FALSE,
  #   linetype = "dashed",
  #   linewidth = 0.5,
  #   color = "#bebdbd"
  # ) +
  scale_fill_manual(values = substrate_colors) +
  scale_color_manual(values = substrate_colors) +
  facet_wrap(~Substrate, ncol = 1) +
  labs(
    x = NULL,
    y = "Relative Abundance (%)",
    title = "Relative Abundance of MAGs across Substrates"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.title.x = element_blank(),
    legend.position = "right",
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
  )

#----------------------------
# 4) Bottom annotation plot for Order labels
#----------------------------
p_order <- ggplot() +
  geom_text(
    data = order_centers,
    aes(x = xmid, y = 1, label = Order_clean),
    size = 3.2, angle = 90
  ) +
  geom_vline(
    data = order_boundaries %>% filter(boundary < max(plot_table2$x_num) + 0.5),
    aes(xintercept = boundary),
    linetype = "dashed",
    linewidth = 0.5,
    color = "#bebdbd"
  ) +
  scale_x_continuous(
    limits = c(0.5, max(plot_table2$x_num) + 0.5),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5)
  )

#----------------------------
# 5) Combine plots
#----------------------------
final_plot <- p_main / p_order +
  plot_layout(heights = c(12, 1.8))

final_plot

ggsave(
  plot = final_plot,
  filename = file.path(output_path, "MAG_singlem_abundance_with_taxid_order.png"),
  width = 15,
  height = 10,
  dpi = 300
)

# ===========================
# Plot assembly metrics
# ===========================

library(scales)
library(forcats)
library(grid)

#---------------------------------
# 1) Prepare x-axis order
#---------------------------------
plot_table3 <- plot_table %>%
  mutate(
    TaxID_clean = 
      sub(" .*$", "", TaxID), #%>%              # remove text after first space
      # gsub("^[a-z]__", "", .),        # remove prefix like s__, g__, etc.
    Order_clean = gsub("^[a-z]__", "", Order),
    x_label = paste0("(", TaxID_clean, ") ", MAG_ID)
    # or use just MAG_ID if cleaner:
    # x_label = MAG_ID
  ) %>%
  arrange(Order_clean, MAG_ID) %>%
  mutate(
    x_label = factor(x_label, levels = unique(x_label))
  )

#---------------------------------
# 2) Make long table for metrics
#---------------------------------
plot_metrics <- plot_table3 %>%
  distinct(
    MAG_ID, x_label, Order_clean,
    completeness, contamination, gc, n50_scaffolds
  ) %>%
  pivot_longer(
    cols = c(completeness, contamination, gc, n50_scaffolds),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    value = ifelse(metric == "n50_scaffolds", value / 1e6, value),
    metric = factor(
      metric,
      levels = c("completeness", "contamination", "gc", "n50_scaffolds"),
      labels = c("Completeness (%)", "Contamination (%)", "GC (%)", "N50 (Mb)")
    )
  )

#---------------------------------
# 3) Order boundaries for separators
#---------------------------------
order_df <- plot_table3 %>%
  distinct(x_label, Order_clean) %>%
  mutate(x_num = as.numeric(x_label))

order_centers <- order_df %>%
  group_by(Order_clean) %>%
  summarise(
    xmin = min(x_num),
    xmax = max(x_num),
    xmid = mean(range(x_num)),
    .groups = "drop"
  )

boundary_df <- order_centers %>%
  mutate(boundary = xmax + 0.5) %>%
  filter(boundary < max(order_df$x_num) + 0.5)

#---------------------------------
# 4) Plot metrics with x-axis on top
#---------------------------------
p_metrics <- ggplot(plot_metrics, aes(x = x_label, y = value)) +
  geom_point(
    aes(color = metric, shape = metric),
    size = 2.8,
    alpha = 0.9,
    stroke = 0.4
  ) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  scale_x_discrete(
    position = "top",
    expand = c(0, 0)
  ) +
  scale_color_manual(
    values = c(
      "Completeness (%)"   = "#1b9e77",
      "Contamination (%)"  = "#d95f02",
      "GC (%)"             = "#7570b3",
      "N50 (Mb)"           = "#e7298a"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Completeness (%)"   = 16,
      "Contamination (%)"  = 17,
      "GC (%)"             = 15,
      "N50 (Mb)"           = 18
    )
  ) +
  labs(
    x = NULL,
    y = NULL#,
    #title = "MAG quality and genome statistics"
  ) +
  theme_bw() +
  theme(
    # axis.text.x.top = element_text(
    #   angle = 90,
    #   hjust = 0,
    #   vjust = 0.5,
    #   size = 7
    # ),
    axis.text.x.top = element_blank(),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 15),
    panel.spacing = unit(0.6, "lines")
  )

final_metrics_plot <-  p_order / p_metrics +
  plot_layout(heights = c(2.2, 4))
final_metrics_plot

# ggsave(
#   plot = final_metrics_plot,
#   filename = file.path(output_path, "MAG_quality_metrics.png"),
#   width = 15,
#   height = 10,
#   dpi = 300
# )


final_combined <- p_main / p_order / p_metrics +
  plot_layout(heights = c(8, 1.8, 6))
final_combined

ggsave(
  plot = final_combined,
  filename = file.path(output_path, "MAG_abundance_and_metrics.png"),
  width = 15,
  height = 16,
  dpi = 300
)
