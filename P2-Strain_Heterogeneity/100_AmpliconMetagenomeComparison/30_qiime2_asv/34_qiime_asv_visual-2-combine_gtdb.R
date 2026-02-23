# Processing RDS data from the combined ASV table
library(tidyverse)
library(janitor)
library(scales)
library(microeco)

# Setting up Work Directory (For non-combined github repo)
#setwd("/home/glbcabria/Workbench/")

# DECLARING AND LOADING INPUTS
inputRDS <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv/dada2_ASV_gtdb_meco.2024.2025.rds"
outputdir <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv/"
dada2_meco <- readRDS(inputRDS)

# Checking Count Distribution
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

# Rarefying the dataset

## Set minimum depth for rarefying based on the histogram and lowest read counts
min_depth <- 39000
set.seed(123)  # for reproducibility
meco_rarefied_list <- lapply(1:3, function(i) {
  tmp <- clone(dada2_meco)
  tmp$rarefy_samples(sample.size = rare_depth)
  return(tmp)
})

# relative abundances + stacked bar at Phylum, top 10
abund_plots <- lapply(seq_along(meco_rarefied_list), function(i) {
  obj <- meco_rarefied_list[[i]]

  # (optional) ensure relative abundance is used (depends on your object state)
  obj$cal_abund()

  p <- obj$plot_bar(
    taxrank = "Phylum",
    topn = 10
  ) +
    ggtitle(paste0("Replicate ", i, " — Phylum (Top 10)"))

  p
})

# Print all plots
for (p in abund_plots) print(p)