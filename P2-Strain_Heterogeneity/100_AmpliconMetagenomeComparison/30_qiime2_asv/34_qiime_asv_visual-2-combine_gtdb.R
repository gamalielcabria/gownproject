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

###############################
# Checking Count Distribution #
###############################

sample_depth <- colSums(dada2_meco$otu_table)
df <- data.frame(   Sample = names(sample_depth),
                    Reads = sample_depth )
lowest10 <- df %>% arrange(Reads) %>%
  slice(1:10)
highest10 <- df %>% arrange(desc(Reads)) %>%
  slice(1:10)
lowest10
highest10

# write.csv(lowest10, paste0(outputdir, "34_qiime_asv_visual-2-combine_gtdb_lowest10.csv"), row.names = FALSE)
# write.csv(highest10, paste0(outputdir, "34_qiime_asv_visual-2-combine_gtdb_highest10.csv"), row.names = FALSE)

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

# ggsave(plot_hist, filename = paste0(outputdir, "34_qiime_asv_visual-2-combine_gtdb_hist.png"), width = 8, height = 6)
# ggsave(plot_hist_lower_half, filename = paste0(outputdir, "34_qiime_asv_visual-2-combine_gtdb_hist_lower_half.png"), width = 8, height = 6)

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

### Replicate 1
# obj <- meco_rarefied_list[[1]]
# obj$cal_abund()
# ta <- trans_abund$new(dataset = obj, taxrank = "Phylum")
# plot_abund_r1 <- ta$plot_bar(others_color = "grey70", facet = c("Year"), xtext_keep = TRUE, legend_text_italic = TRUE, barwidth = 1) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))

# ### Replicate 2
# obj <- meco_rarefied_list[[1]]
# obj$cal_abund()
# ta <- trans_abund$new(dataset = obj, taxrank = "Phylum")
# plot_abund_r2 <- ta$plot_bar(others_color = "grey70", facet = c("Year"), xtext_keep = TRUE, legend_text_italic = TRUE, barwidth = 1) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))

# ### Replicate 3
# obj <- meco_rarefied_list[[1]]
# obj$cal_abund()
# ta <- trans_abund$new(dataset = obj, taxrank = "Phylum")
# plot_abund_r3 <- ta$plot_bar(others_color = "grey70", facet = c("Year"), xtext_keep = TRUE, legend_text_italic = TRUE, barwidth = 1) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))

# ggsave(plot_abund_r1, filename = paste0(outputdir, "34_qiime_asv_visual-2-combine_gtdb_abund_r1.png"), width = 10, height = 6)
# ggsave(plot_abund_r2, filename = paste0(outputdir, "34_qiime_asv_visual-2-combine_gtdb_abund_r2.png"), width = 10, height = 6)
#ggsave(plot_abund_r3, filename = paste0(outputdir, "34_qiime_asv_visual-2-combine_gtdb_abund_r3.png"), width = 10, height = 6)


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

### --- sanity check ---
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
#palettes <- c("#6ea577", "#D7B18C", "#F0A390", "#AAB0D4", "#D0A5D8", "#ECA8BE", "#BDD676", 
#"#D7E151", "#F28E2B", "#F6D383", "#DFC9AB", "#C2C2C2", "grey70")

base_cols <- colorRampPalette(brewer.pal(8, "Set2"))(length(tax_no_others))
soft_cols <- lighten(base_cols, amount = 0.2)
pal_main <- setNames(soft_cols, tax_no_others)
pal <- c(pal_main, Others = "grey70")
pal["Patescibacteria"] <- "#F28E2B" # Force them to be
pal["Actinomycetota"] <- "#6ea577"

plot_abund_r1 <- ggplot(plot_df, aes(x = Well, y = Abundance, fill = Taxon2)) +
  geom_col(width = 1) +
  facet_wrap(~Year, scales = "free_x") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
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
ggsave(plot_abund_r1, filename = paste0(outputdir, "34_qiime_asv_visual-2-combine_gtdb_abund_r1.png"), width = 10, height = 6)



############################
# ALPHA AND BETA DIVERSITY #
############################

## ALPHA DIVERSITY
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

# combined_alpha <- combined_alpha %>%
#   left_join(meta %>% select(Sample, Well, Year), by = "Sample") #%>%
#   mutate(
#     Well = as.factor(Well),
#     Year = as.factor(Year)
#   )

alpha_metrics <- c("Observed", "Chao1", "Shannon", "Simpson", "Pielou")

alpha_summary <- combined_alpha%>%
  select(Sample, Well, Year, Measure, Value, Replicate) %>%
  filter(Measure %in% alpha_metrics) %>%
  group_by(Well, Year, Measure) %>%
  arrange(Sample) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD   = sd(Value, na.rm = TRUE),
    SE   = SD / sqrt(n()),
    .groups = "drop"
  )
alpha_summary


# Visualize
plot_alpha <- ggplot(alpha_summary, aes(x = Well, y = Mean)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.25) +
  geom_point(size = 2) +
  # show replicate-level variation
  geom_point(
    data = combined_alpha,
    aes(x = Well, y = Value),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.15, height = 0),
    alpha = 0.3,
    size = 1
  ) +
  facet_grid(Measure ~ Year, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    strip.text.y = element_text(face = "bold"),
    panel.spacing = unit(0.4, "lines")
  ) +
  labs(
    title = paste0("Alpha diversity across rarefaction replicates (n = ",
                   length(meco_rarefied_list), ")"),
    x = "Well",
    y = "Mean ± SD"
  )

plot_alpha

# ggsave(plot_alpha, filename = paste0(outputdir, "34_qiime_asv_visual-2-combine_gtdb_alpha.png"), width = 10, height = 6)


## BETA DIVERISTY
obj <- meco_rarefied_list[[1]]
obj$cal_betadiv(measure="bray")
tbd <- trans_beta$new(dataset = obj, group = "Well", measure = "bray")

### PCoA
tbd$cal_ordination(method = "PCoA")
class(tbd$res_ordination)
plot_beta <-tbd$plot_ordination()
plot_beta

ggsave(plot_beta, filename = paste0(outputdir, "34_qiime_asv_visual-2-combine_gtdb_beta-PCoA.png"), width = 10, height = 6)

### NMDS
tbd$cal_ordination(method = "DCA")
class(tbd$res_ordination)
plot_beta <-tbd$plot_ordination()
plot_beta

ggsave(plot_beta, filename = paste0(outputdir, "34_qiime_asv_visual-2-combine_gtdb_beta-DCA.png"), width = 10, height = 6)

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

ggsave(plot_beta, filename = paste0(outputdir, "34_qiime_asv_visual-2-combine_gtdb_beta-NMDS.png"), width = 10, height = 6)
