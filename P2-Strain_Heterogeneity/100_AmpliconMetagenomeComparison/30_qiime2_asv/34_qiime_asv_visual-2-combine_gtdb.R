# Processing RDS data from the combined ASV table
library(tidyverse)
library(janitor)
library(scales)
library(microeco)
library(ggh4x)
library(RColorBrewer)
library(patchwork)

# Setting up Work Directory (For non-combined github repo)
#setwd("/home/glbcabria/Workbench/")

# DECLARING AND LOADING INPUTS
inputRDS <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv/dada2_ASV_gtdb_meco.2024.2025.rds"
outputdir <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv/"
dada2_meco <- readRDS(inputRDS)
output_tag <- tools::file_path_sans_ext(basename(inputRDS))

###############################
# Checking Count Distribution #
###############################

sample_depth <- colSums(dada2_meco$otu_table)
df <- data.frame(   Sample = names(sample_depth),
                    Reads = sample_depth )
lowest10 <- df %>% arrange(Reads) %>%  slice(1:10)
highest10 <- df %>% arrange(desc(Reads)) %>%  slice(1:10)
lowest10
highest10

write.csv(lowest10, paste0(outputdir, "34_qiime_asv_visual-2__", output_tag, "__gtdb_lowest10.csv"), row.names = FALSE)
write.csv(highest10, paste0(outputdir, "34_qiime_asv_visual-2__", output_tag, "__gtdb_highest10.csv"), row.names = FALSE)

## Plot count histogram
cutoff <- quantile(df$Reads, 0.5)
bin_w <- cutoff / 20   # more bins → narrower bars

plot_hist_lower_half <- ggplot(df, aes(x = Reads)) +
  geom_histogram(binwidth = bin_w,
                 fill = "grey70",
                 color = "black") +
  coord_cartesian(xlim = c(0, cutoff)) +
  theme_bw() +
  labs(title = "Lower Tail of Read Depth Distribution",
       x = "Total Reads",
       y = "Number of Samples")

plot_hist <- ggplot(df, aes(x = Reads)) +
  geom_histogram(bins = 70, fill = "grey80", color = "black") +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),   # more tick marks
    labels = label_number(big.mark = ",")               # no scientific notation
  ) +
  theme_bw() +
  labs(title = "Log-Scaled Read Depth Distribution",
       x = "Total Reads",
       y = "Number of Samples") +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold")
  )

ggsave(plot_hist, filename = paste0(outputdir, "34_qiime_asv_visual-2__", output_tag, "__gtdb_hist.png"), width = 8, height = 6)
ggsave(plot_hist_lower_half, filename = paste0(outputdir, "34_qiime_asv_visual-2__", output_tag, "__gtdb_hist_lower_half.png"), width = 8, height = 6)

#########################
# Rarefying the dataset #
#########################

## Set minimum depth for rarefying based on the histogram and lowest read counts
rare_depth <- 25000
set.seed(913)  # for reproducibility
meco_rarefied_list <- lapply(1:3, function(i) {
  tmp <- clone(dada2_meco)
  tmp$rarefy_samples(sample.size = rare_depth)
  return(tmp)
})

#########################################
# Abundance Visualization in Replicates #
#########################################

# Abundance plot V2
taxrank_use <- "Phylum"
topn_use <- 12

### Change this 
obj <- meco_rarefied_list[[1]]
obj$cal_abund()

ta <- trans_abund$new(dataset = obj, taxrank = taxrank_use)

### From Microeco to matrix
abund_df <- as.data.frame(ta$data_abund, stringsAsFactors = FALSE) %>% as_tibble()
abund_df <- abund_df %>%
  mutate(
    Abundance = as.numeric(Abundance),
    Year = as.factor(Year),
    Well = as.factor(Well),
    Taxonomy = as.character(Taxonomy)
  )
meta <- obj$sample_table %>% as.data.frame() %>%
  tibble::rownames_to_column("Sample")

### ~~ SANITY CHECK ~~
stopifnot(all(abund_df$Sample %in% meta$Sample))
stopifnot(all(c("Well","Year") %in% colnames(meta)))

### Pick top taxa across all samples by mean abundance
top_taxa <- abund_df %>%
  group_by(Taxonomy) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = topn_use) %>%
  pull(Taxonomy)

# Collapse to Others, then sum within Well-Year (so one bar per Well)
plot_df <- abund_df %>%
  mutate(Taxon2 = if_else(Taxonomy %in% top_taxa, Taxonomy, "Others")) %>%
  group_by(Year, Well, Taxon2, Sample) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(Taxon2 = factor(Taxon2, levels = c(sort(setdiff(unique(Taxon2), "Others")), "Others"))) # Force "Others" to be last in legend

# Make palette with Others = grey70
#palettes <- c("#6ea577", "#D7B18C", "#F0A390", "#AAB0D4", "#D0A5D8", "#ECA8BE", "#BDD676", "#D7E151", "#F28E2B", "#F6D383", "#DFC9AB", "#C2C2C2", "grey70")

base_cols <- colorRampPalette(brewer.pal(8, "Set2"))(length(tax_no_others))
soft_cols <- lighten(base_cols, amount = 0.2)
pal_main <- setNames(soft_cols, tax_no_others)
pal <- c(pal_main, Others = "grey70")
pal["Patescibacteria"] <- "#F28E2B" # Force them to be
pal["Actinomycetota"] <- "#6ea577"

### Removing Replicates
pick_suffixes <- c("a")
unique_df <- plot_df %>%
  mutate(   BaseSample = sub("([A-Za-z]+\\d+).*", "\\1", Sample),
            Rep = sub("^[A-Za-z]+\\d+", "", Sample)) %>%
  filter(Rep %in% pick_suffixes | Rep == "") %>%
  select(-Rep)

plot_abund_r1 <- ggplot(unique_df, aes(x = Well, y = Abundance, fill = Taxon2)) +
  geom_col(width = 1) +
  # facet_wrap(~Year, ncol = 1)+#, scales = "free_x") +
  facet_grid(rows = vars(Year), scales = "free_x", switch = "y") +
  scale_y_continuous(labels = percent_format(accuracy = 1, scale = 1)) +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.text = element_text(face = "italic"),
    panel.spacing = unit(0.4, "lines")
  ) +
  labs(
    title = paste0("Taxonomic composition (", taxrank_use, ") — Replicate 1"),
    x = "Well", y = "Relative abundance",
    fill = taxrank_use
  )

plot_abund_r1
ggsave(plot_abund_r1, filename = paste0(outputdir, "34_qiime_asv_visual-2__", output_tag, "__gtdb_abund_r1.png"), width = 10, height = 6)


# COMMON SAMPLES ONLY

## 1) Find wells common to 2024 and 2025
years_interest <- c("2024", "2025")   # make sure Year in unique_df is a factor with these levels

wells_common <- unique_df %>%
  filter(Year %in% years_interest) %>%
  distinct(Year, Well) %>%
  count(Well, name = "n_years") %>%
  filter(n_years == length(years_interest)) %>%
  pull(Well)
message("Common wells (", length(wells_common), "): ", paste(sort(as.character(wells_common)), collapse = ", "))

## 2) Filter data to those wells and years; fix factor orders for clean facets
unique_df_common <- unique_df %>%
  filter(Year %in% years_interest, Well %in% wells_common) %>%
  mutate(
    Year = factor(Year, levels = years_interest),
    # consistent Well order across facets; change to a custom order if you prefer
    Well = fct_relevel(Well, sort(levels(factor(wells_common))))
  ) %>%
  droplevels()

## 3) Plot: only wells common to both years (facet labels on the left)
plot_abund_common <- ggplot(unique_df_common, aes(x = Well, y = Abundance, fill = Taxon2)) +
  geom_col(width = 1) +
  facet_grid(rows = vars(Year), scales = "fixed", switch = "y") +  # fixed keeps the same wells/order across years
  scale_y_continuous(labels = percent_format(accuracy = 1, scale = 1)) +
  scale_fill_manual(values = pal) +
  theme_bw() +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", colour = NA),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.text = element_text(face = "italic"),
    panel.spacing = unit(0.4, "lines")
  ) +
  labs(
    title = paste0("Taxonomic composition (", taxrank_use, ") — Wells common to 2024 & 2025 (Replicate 1)"),
    x = "Well", y = "Relative abundance",
    fill = taxrank_use
  )
plot_abund_common

ggsave(plot_abund_common,
       filename = paste0(outputdir, "34_qiime_asv_visual-2__", output_tag, "__gtdb_abund_common.png"),
       width = 6, height = 6, dpi = 300)


############################
# ALPHA AND BETA DIVERSITY #
############################
##-----------------
## ALPHA DIVERSITY
##-----------------
combined_alpha <- purrr::map_dfr(seq_along(meco_rarefied_list), function(i) {
  obj <- meco_rarefied_list[[i]]
  obj$cal_alphadiv()
  tad <- trans_alpha$new(dataset = obj)
  tad$data_alpha %>%
    as.data.frame() %>%
    # mutate(ItemID = paste0("Rep", i, "_", row_number())) %>%
    #rownames_to_column("Sample") %>%   # keeps sample IDs explicit
    mutate(Replicate = i)
})

# -- sanity check --
combined_alpha %>% count(Replicate)
combined_alpha %>% glimpse()
alpha_metrics <- c("Observed", "Shannon", "Simpson", "Pielou")

combined_alpha <- combined_alpha %>% mutate( Well = as.factor(Well), Year = as.factor(Year) )
combined_alpha_filtered <- combined_alpha %>% filter(Measure %in% alpha_metrics)

### Set which metrics to visualize
alpha_summary <- combined_alpha %>%
  select(Sample, Well, Year, Measure, Value, Replicate) %>%
  filter(Measure %in% alpha_metrics) %>%
  group_by(Well, Year, Measure) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD   = sd(Value, na.rm = TRUE),
    N    = dplyr::n(),
    SE   = SD / sqrt(N),
    .groups = "drop"
  )

# Reusable plotting function (no custom attributes)
plot_alpha_metric <- function(metric_name, alpha_summary, combined_alpha_filtered, n_reps) {
  ggplot(
    data = alpha_summary %>% filter(Measure == metric_name),
    aes(x = Well, y = Mean)
  ) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.25) +
    geom_point(size = 2) +
    # replicate-level points
    geom_point(
      data = combined_alpha_filtered %>% filter(Measure == metric_name),
      aes(x = Well, y = Value),
      inherit.aes = FALSE,
      position = position_jitter(width = 0.15, height = 0),
      alpha = 0.3, size = 1
    ) +
    # Year stacked vertically (on top of each other)
    facet_grid(rows = vars(Year)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      strip.text.y = element_text(face = "bold"),
      panel.spacing = unit(0.4, "lines")
    ) +
    labs(
      title = paste0(metric_name, " (n = ", n_reps, ")"),
      x = "Well", y = "Mean ± SD"
    )
}

# Save each metric separately
n_reps <- length(meco_rarefied_list)

p_metric <- plot_alpha_metric(
    metric_name = m,
    alpha_summary = alpha_summary, combined_alpha_filtered = combined_alpha_filtered, n_reps = n_reps
  )



plots_by_metric <- map(
  alpha_metrics,
  ~ plot_alpha_metric(
      metric_name = .x,
      alpha_summary = alpha_summary,
      combined_alpha_filtered = combined_alpha_filtered,
      n_reps = n_reps
    )
) %>%
  set_names(alpha_metrics)

combined_patchwork <- wrap_plots(plots_by_metric, ncol = length(alpha_metrics)) 

ggsave(
  filename = file.path(outputdir, paste0("34_qiime_asv_visual-2__", output_tag, "__gtdb_alpha_combined.png")),
  plot = combined_patchwork,
  width = 12, height = 6, dpi = 300
)

###-------------------------------------------
### COMBINED ALPHA DIVERSITY ACROSS THE YEARS
###-------------------------------------------
# Ensure types are consistent
combined_alpha <- combined_alpha %>% mutate(Well = as.factor(Well), Year = as.factor(Year))

# 1) Wells that occur in BOTH 2024 and 2025
common_wells <- combined_alpha %>%
  distinct(Well, Year) %>%
  mutate(Year = as.character(Year)) %>%
  filter(Year %in% c("2024", "2025")) %>%
  pivot_wider(
    names_from = Year, values_from = Year,
    values_fill = NA, values_fn = length
  ) %>%
  filter(!is.na(`2024`), !is.na(`2025`)) %>%
  pull(Well)
n_common_wells <- length(common_wells)
message("Wells present in both 2024 and 2025: ", n_common_wells)

# 2) Filter replicate-level alpha for only those common wells + selected metrics
combined_alpha_bothyrs <- combined_alpha %>%
  filter(Well %in% common_wells,
         Measure %in% alpha_metrics) %>%
  droplevels()

# 3) Recompute summary on the intersection set
alpha_summary_bothyrs <- combined_alpha_bothyrs %>%
  select(Sample, Well, Year, Measure, Value, Replicate) %>%
  group_by(Well, Year, Measure) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD   = sd(Value, na.rm = TRUE),
    N    = dplyr::n(),
    SE   = SD / sqrt(N),
    .groups = "drop"
  )

# 4) Plot helper (same style as yours; adds subtitle + wells count)
plot_alpha_metric_wells_bothyrs <- function(metric_name, alpha_summary, combined_alpha_filtered, n_reps, n_common_wells) {
  ggplot(
    data = alpha_summary %>% filter(Measure == metric_name),
    aes(x = Well, y = Mean)
  ) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.25) +
    geom_point(size = 2) +
    # replicate-level points (jittered)
    geom_point(
      data = combined_alpha_filtered %>% filter(Measure == metric_name),
      aes(x = Well, y = Value),
      inherit.aes = FALSE,
      position = position_jitter(width = 0.15, height = 0),
      alpha = 0.3, size = 1
    ) +
    facet_grid(rows = vars(Year)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      strip.text.y = element_text(face = "bold"),
      panel.spacing = unit(0.4, "lines")
    ) +
    labs(
      title = paste0(metric_name, " (n reps = ", n_reps, ")"),
      x = "Well", y = "Mean ± SD"
    )
}

# 5) Build patchwork of all metrics for the intersection set
n_reps <- length(meco_rarefied_list)

plots_by_metric_wells_bothyrs <- map(
  alpha_metrics,
  ~ plot_alpha_metric_wells_bothyrs(
      metric_name = .x,
      alpha_summary = alpha_summary_bothyrs,
      combined_alpha_filtered = combined_alpha_bothyrs,
      n_reps = n_reps,
      n_common_wells = n_common_wells
    )
) %>% set_names(alpha_metrics)

combined_patchwork_wells_bothyrs <- wrap_plots(plots_by_metric_wells_bothyrs, ncol = length(alpha_metrics)) +  
plot_annotation(
    title = "Alpha Diversity Metrics (Common Wells in 2024 & 2025)",
    subtitle = paste0("Wells in both years included (n = ", n_common_wells, ")" ),
    caption = "Error bars = SD across rarefaction replicates"
  )
ggsave(
  filename = file.path(outputdir, paste0("34_qiime_asv_visual-2__", output_tag, "__gtdb_alpha_commonWells_2024_2025.png")),
  plot = combined_patchwork_wells_bothyrs,
  width = 12, height = 6, dpi = 300
)



##----------------
## BETA DIVERISTY
##----------------
obj <- meco_rarefied_list[[1]]
obj$cal_betadiv(measure="bray")
obj$sample_table$Year <- as.factor(obj$sample_table$Year)
tbd <- trans_beta$new(dataset = obj, group = "Well", measure = "bray")

### PCoA
tbd$cal_ordination(method = "PCoA")
class(tbd$res_ordination)
plot_beta <-tbd$plot_ordination(plot_color="Year")
plot_beta

ggsave(plot_beta, filename = paste0(outputdir, "34_qiime_asv_visual-2__", output_tag, "__gtdb_beta-PCoA.png"), width = 10, height = 6)

### NMDS
tbd$cal_ordination(method = "DCA")
class(tbd$res_ordination)
plot_beta <-tbd$plot_ordination(plot_color="Year")
plot_beta

ggsave(plot_beta, filename = paste0(outputdir, "34_qiime_asv_visual-2__", output_tag, "__gtdb_beta-DCA.png"), width = 10, height = 6)

### NMDS
tbd$cal_ordination(method = "NMDS", ncomp = 3)
class(tbd$res_ordination)
nmds_df <- as.data.frame(tbd$res_ordination$scores)
nmds_df$Sample <- rownames(nmds_df)
stress_val <- tbd$res_ordination$model$stress
plot_df <- nmds_df

p12 <- ggplot(plot_df, aes(x = MDS1, y = MDS2, color = Year)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "NMDS (k = 3): Axes 1 vs 2")
p13 <- ggplot(plot_df, aes(x = MDS1, y = MDS3, color = Year)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "NMDS (k = 3): Axes 1 vs 3")
p23 <- ggplot(plot_df, aes(x = MDS2, y = MDS3, color = Year)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "NMDS (k = 3): Axes 2 vs 3")

plot_beta <- (p12 | p13) / (p23 | plot_spacer()) +
  plot_annotation(
    caption = paste0("NMDS stress = ", round(stress_val, 4))
  ) &
  theme(
    plot.caption = element_text(size = 10, face = "bold"),
    plot.caption.position = "plot"  # keeps it at the bottom edge of the full patchwork
  )

ggsave(plot_beta, filename = paste0(outputdir, "34_qiime_asv_visual-2__", output_tag, "__gtdb_beta-NMDS.png"), width = 10, height = 6)
