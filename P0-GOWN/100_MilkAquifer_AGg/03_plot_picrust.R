#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

infile <- "03_KO_process_by_sample_taxon_fallback_summary.tsv"

outdir <- "04_plots_KO_taxon_dotplot"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Processes to plot
# Set to NULL to plot all
# -----------------------------
process_keep <- c(
  "denitrification",
  #"DNRA",
  "methanotrophy",
  "methanogenesis",
  #"nitrification",
  #"anammox",
  "sulfate_reduction",
  "sulfur_oxidation",
  "iron_cycling"
)

# -----------------------------
# Sample renaming
# -----------------------------
sample_map <- c(
  "DMO-GOWN-bact-2023-G21-20230321" = "GOWN 211 - Smith Coulee",
  "DMO-GOWN-bact-2023-G22-20230321" = "GOWN 213 - Milk River 85-2",
  "DMO-GOWN-bact-2023-G23-20230321" = "GOWN 212 - Pine Coulee"
)

# Samples to plot
# Set to NULL to plot all samples
sample_keep <- c(
  "GOWN 211 - Smith Coulee",
  "GOWN 213 - Milk River 85-2"
)

# -----------------------------
# Read and filter
# -----------------------------
df <- read_tsv(infile, show_col_types = FALSE) %>%
  filter(!is.na(Process), !is.na(Taxon_label)) %>%
  mutate(
    sample_label = recode(sample, !!!sample_map)
  )

if (!is.null(process_keep)) {
  df <- df %>%
    filter(Process %in% process_keep)
}

if (!is.null(sample_keep)) {
  df <- df %>%
    filter(sample_label %in% sample_keep)
}

# -----------------------------
# Step 1: sum ASVs within same Gene
# -----------------------------
df_gene <- df %>%
  group_by(Process, sample_label, Taxon_rank, Taxon_label, Gene) %>%
  summarise(
    gene_taxon_rel_abun = sum(total_taxon_rel_abun, na.rm = TRUE),
    gene_norm_contribution = sum(total_norm_contribution, na.rm = TRUE),
    n_ASVs = sum(n_ASVs, na.rm = TRUE),
    KOs = paste(sort(unique(KO)), collapse = "; "),
    .groups = "drop"
  )

# -----------------------------
# Step 2: mean across genes within process
# -----------------------------
df_plot <- df_gene %>%
  group_by(Process, sample_label, Taxon_rank, Taxon_label) %>%
  summarise(
    plot_rel_abun = mean(gene_taxon_rel_abun, na.rm = TRUE),
    plot_norm_contribution = mean(gene_norm_contribution, na.rm = TRUE),
    n_genes = n_distinct(Gene),
    genes = paste(sort(unique(Gene)), collapse = "; "),
    KOs = paste(sort(unique(KOs)), collapse = "; "),
    n_ASVs = sum(n_ASVs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Taxon_label_plot = paste0(Taxon_label, " [", Taxon_rank, "]")
  )

min_rel_abun <- 1

df_plot <- df_plot %>%
  filter(plot_rel_abun > min_rel_abun)

write_tsv(df_plot, file.path(outdir, "00_KO_process_taxon_dotplot_plotting_table.tsv"))

sample_order <- sample_keep

df_plot <- df_plot %>%
  mutate(
    sample_label = factor(sample_label, levels = sample_order),
    Taxon_label_plot = fct_reorder(Taxon_label_plot, plot_rel_abun, .fun = max)
  )

pastel_cols <- c(
  "#F9A5A4", "#92C5DE", "#BEAED4", "#FDB462", "#ffff89",
  "#B3DE69", "#80B1D3", "#FCCDE5", "#BC80BD", "#8DD3C7",
  "#CCEBC5", "#FFED6F", "#BEBADA", "#FB8072", "#D9D9D9"
)

# -----------------------------
# Combined plot: process facets on TOP
# -----------------------------
p <- ggplot(df_plot, aes(x = sample_label, y = Taxon_label_plot)) +
  geom_point(aes(size = plot_rel_abun, color = Process), alpha = 0.85) +
  facet_wrap(~ Process, scales = "fixed", nrow = 1) +
  scale_color_manual(values = rep(pastel_cols, length.out = n_distinct(df_plot$Process))) +
  scale_size_continuous(name = "Relative abundance") +
  labs(
    x = "Sample",
    y = "Taxon label",
    color = "Function / process",
    title = "PICRUSt2 KO-linked taxon relative abundance by function",
    subtitle = "ASVs summed within each gene; multiple genes averaged within each process"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 12),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(outdir, "01_KO_process_taxon_dotplot_all_functions.png"),
  p,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

message("Done.")
message("Plot written to: ", outdir)