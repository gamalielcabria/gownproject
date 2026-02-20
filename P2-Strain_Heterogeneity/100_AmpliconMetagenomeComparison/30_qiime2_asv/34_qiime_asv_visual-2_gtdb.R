# Processing RDS data from the combined ASV table
library(tidyverse)
library(janitor)
library(microeco)

# Setting up Work Directory (For non-combined github repo)
setwd("/home/glbcabria/Workbench/")

# DECLARING AND LOADING INPUTS
inputRDS <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv/dada2_ASV_gtdb_meco.2024.2025.rds"

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

#
cutoff <- quantile(df$Reads, 0.25)

ggplot(df, aes(x = Reads)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black") +
  coord_cartesian(xlim = c(0, cutoff)) +
  theme_bw() +
  labs(title = "Lower Tail of Read Depth Distribution",
       x = "Total Reads",
       y = "Number of Samples")
