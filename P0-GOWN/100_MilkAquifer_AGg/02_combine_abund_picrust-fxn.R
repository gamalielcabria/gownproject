#!/usr/bin/env Rscript

library(tidyverse)

# Inputs
ko_contrib_file <- "picrust2_out/KO_contrib_out/pred_metagenome_contrib.tsv.gz"
ko_map_file     <- "ko_process_mapping.tsv"
taxa_file       <- "taxa.csv"
fasta_file      <- "picrust_ASV.fasta"

# Optional: specific samples to keep
samples_keep <- c(
  "DMO-GOWN-bact-2023-G21-20230321",
  "DMO-GOWN-bact-2023-G22-20230321",
  "DMO-GOWN-bact-2023-G23-20230321"
)

# -----------------------------
# 1. Read ASV fasta map
# -----------------------------
fasta_lines <- readLines(fasta_file)

asv_ids <- fasta_lines[grepl("^>", fasta_lines)] %>%
  sub("^>", "", .)

seqs <- fasta_lines[!grepl("^>", fasta_lines)]

asv_map <- tibble(
  ASV = asv_ids,
  sequence = seqs
)

# -----------------------------
# 2. Read taxonomy
# -----------------------------
taxa <- read_csv(taxa_file, show_col_types = FALSE) %>%
  rename(sequence = 1)

# -----------------------------
# 3. Read KO process mapping
# -----------------------------
ko_map <- read_tsv(ko_map_file, show_col_types = FALSE) %>%
  filter(!is.na(KO), KO != "")

# -----------------------------
# 4. Read PICRUSt2 KO contribution table
# -----------------------------
ko_contrib <- read_tsv(ko_contrib_file, show_col_types = FALSE) %>%
  mutate(
    KO = sub("^ko:", "", `function`),
    ASV = taxon
  ) %>%
  filter(sample %in% samples_keep)

#-----------------------------
# 5. Join everything
# -----------------------------
combined <- ko_contrib %>%
  left_join(ko_map, by = "KO") %>%
  left_join(asv_map, by = "ASV") %>%
  left_join(taxa, by = "sequence") %>%
  select(
    sample,
    ASV,
    KO,
    Process,
    Gene,
    Description,
    Confidence,
    taxon_abun,
    taxon_rel_abun,
    genome_function_count,
    taxon_function_abun,
    taxon_rel_function_abun,
    norm_taxon_function_contrib,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    sequence
  )

# Save full joined table
write_tsv(combined, "01_KO_ASV_taxa_process_combined.tsv")

# Save only mapped target functions
combined_mapped <- combined %>%
  filter(!is.na(Process))

write_tsv(combined_mapped, "02_KO_ASV_taxa_process_TARGETS_ONLY.tsv")

# Summarize by sample, process, genus
summary_genus <- combined_mapped %>%
  group_by(sample, Process, Gene, Genus) %>%
  summarise(
    total_taxon_abun = sum(taxon_abun, na.rm = TRUE),
    total_taxon_function_abun = sum(taxon_function_abun, na.rm = TRUE),
    total_norm_contribution = sum(norm_taxon_function_contrib, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(sample, Process, desc(total_norm_contribution))

write_tsv(summary_genus, "03_KO_process_by_sample_genus_summary.tsv")
