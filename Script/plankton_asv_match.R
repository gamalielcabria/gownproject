# R script for determining what taxa is associated with matched query sequences
# Usage: Rscript plankton_asv_match.R <path_to_taxonomy_file> <path_to_mmseqs_output_file>

library(tidyverse)
# Input files for testing
input <- "/home/gam/workbench/GMPollard/OldConsensusSeqsForTest2025Dec08_v_SOZoops18SOperon-ReferenceSeqs.output.m8"  # Path to mmseqs output file 
taxon <- "/home/gam/workbench/GMPollard/taxonomy_results.csv" #Path to where the csv converted taxonomy file is located

# Taking in inputs for SLURM
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) != 2) {
#   stop("Usage: Rscript annotate_alignment.R taxonomy_results.csv alignment.m8\n", call. = FALSE)
# }
# taxon <- args[1]
# input <- args[2]

# Load required Libraries 
suppressPackageStartupMessages({
  library(tidyverse)
})

# Open taxonomy file and input alignment file
taxon_df <- read_csv(taxon, skip_empty_rows = TRUE, col_names = TRUE, show_col_types = FALSE) 
alignment_df <- read_tsv(
  input,
  col_names = c("query","target","pident","alnlen","mismatch","gapopen",
    "qstart","qend","tstart","tend","evalue","bitscore"),
  show_col_types = FALSE
)

# Process the taxonomy df to have filled columns and select relevant columns
tax_prefix <- c( #Kingdom = "k__",
Phylum = "p__",Subphylum = "sp__",Infraphylum = "ip__",Parvphylum = "pp__",Gigaclass = "gc__",Superclass = "sc__",Class = "c__",
Subclass = "sb__",Infraclass = "ic__",Subterclass = "stc__",Superorder = "so__",Order = "o__",Family = "f__",Genus = "g__",Species = "s__"
)

taxon_df_filled <- taxon_df %>%
  select(Name, all_of(names(tax_prefix))) %>%
  mutate(across(-Name, ~ ifelse(is.na(.), paste0(tax_prefix[cur_column()], "unclassified"), 
    paste0(tax_prefix[cur_column()], .))))

## Fill down the taxonomy levels
tax_cols <- names(tax_prefix)
strip_prefix <- function(x) sub("^[^_]+__", "", x)

fill_from_parent <- function(row) {
  core <- strip_prefix(row)
  is_uncl <- tolower(core) == "unclassified"

  for (i in which(is_uncl)) {
    if (i == 1) next
    j <- max(which(
      !is.na(core[1:(i-1)]) &
      core[1:(i-1)] != "" &
      tolower(core[1:(i-1)]) != "unclassified"
    ), na.rm = TRUE)

    if (is.finite(j)) row[i] <- row[j]
  }
  row
}

taxon_df_filled2 <- taxon_df_filled
taxon_df_filled2[tax_cols] <- t(apply(taxon_df_filled2[tax_cols], 1, fill_from_parent))

# Combine alignment with taxonomy information
combined_df <- alignment_df %>%
  left_join(taxon_df_filled2, by = c("target" = "Name"))

# Taxonomy columns (use only those that exist)
tax_cols_all <- c("Domain","Kingdom","Phylum","Subphylum","Infraphylum","Parvphylum","Gigaclass","Superclass",
                  "Class","Subclass","Infraclass","Subterclass","Superorder","Order",
                  "Family","Genus","Species")
tax_cols <- intersect(tax_cols_all, names(combined_df))

# Shared taxonomy across ALL hits for a query (only ranks where all hits agree)
# get LCA function
get_lca <- function(df_q, tax_cols) {
  if (length(tax_cols) == 0) return(tibble(lca_rank = NA_character_, lca_taxon = NA_character_))
  agreed <- map_lgl(tax_cols, function(col) {
    vals <- df_q[[col]]
    vals <- vals[!is.na(vals) & vals != ""]
    length(vals) == nrow(df_q) && dplyr::n_distinct(vals) == 1
  })
  if (!any(agreed)) return(tibble(lca_rank = NA_character_, lca_taxon = NA_character_))
  idx <- max(which(agreed))
  rank <- tax_cols[[idx]]
  taxon <- unique(df_q[[rank]])[1]
  tibble(lca_rank = rank, lca_taxon = taxon)
}

# 1) Summary per query (safe summarise)
per_query_summary <- combined_df %>%
  group_by(query) %>%
  summarise(
    n_hits = n(),
    avg_pident = mean(pident, na.rm = TRUE),
    .groups = "drop"
  )

# 2) best hit per query (bitscore desc, then pident desc) 
# (replace NA bitscore/pident with -Inf so they sort to the bottom)
best_hit <- combined_df %>%
  mutate(
    bitscore_rank = if_else(is.na(bitscore), -Inf, bitscore),
    pident_rank   = if_else(is.na(pident),   -Inf, pident)
  ) %>%
  arrange(query, desc(bitscore_rank), desc(pident_rank)) %>%
  group_by(query) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    query,
    best_target   = target,
    best_bitscore = bitscore,
    best_pident   = pident,
    best_alnlen   = alnlen,
    best_mismatch = mismatch
  )

#  3) LCA per query 
lca_tbl <- combined_df %>%
  group_by(query) %>%
  group_modify(~ get_lca(.x, tax_cols)) %>%
  ungroup()

# Final table (one row per query)
per_query_best_lca <- per_query_summary %>%
  left_join(best_hit, by = "query") %>%
  left_join(lca_tbl,  by = "query") %>%
  arrange(desc(n_hits), desc(avg_pident))

write_csv(per_query_best_lca, file = "plankton_asv_mmseqs_taxonomy_summary.csv")
