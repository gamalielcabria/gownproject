#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

# =========================
# Paths
# =========================
base_dir <- "/work/ebg_lab/eb/AIWIP_GOWN_amplicons/400_TimeSeries/20_abund_analysis"

# =========================
# Helpers
# =========================

# Convert counts to relative abundance by sample
# Assumes sample columns start after metadata columns
make_relabund <- function(df, meta_cols) {
  out <- copy(df)

  sample_cols <- setdiff(colnames(out), meta_cols)

  for (s in sample_cols) {
    total <- sum(out[[s]], na.rm = TRUE)
    if (is.na(total) || total == 0) {
      out[[s]] <- 0
    } else {
      out[[s]] <- out[[s]] / total
    }
  }

  return(out)
}

# Fill the taxa based on the parent
fill_from_parent <- function(row) {
  core <- str_remove(row, "^[a-z]__")
  is_uncl <- tolower(core) == "unclassified"

  for (i in which(is_uncl)) {
    if (i == 1) next

    prev_core <- str_remove(row[1:(i-1)], "^[a-z]__")
    prev_idx <- which(
      !is.na(row[1:(i-1)]) &
      row[1:(i-1)] != "" &
      tolower(prev_core) != "unclassified"
    )

    if (length(prev_idx) > 0) {
      j <- max(prev_idx)
      row[i] <- row[j]
    }
  }
  row
}

# Read and transpose otu table
# Input format:
# - first row = "" then sequences
# - first column = sample names
# - remaining entries = counts
read_transpose_otu <- function(file) {
  otu_raw <- fread(file, check.names = FALSE)

  first_col <- colnames(otu_raw)[1]
  setnames(otu_raw, first_col, "Sample")

  otu_long <- melt(
    otu_raw,
    id.vars = "Sample",
    variable.name = "Sequence",
    value.name = "Count"
  )

  otu_long[, Count := as.numeric(Count)]

  otu_wide <- dcast(
    otu_long,
    Sequence ~ Sample,
    value.var = "Count",
    fill = 0
  )

  return(otu_wide)
}

# =========================
# Find run folders
# =========================
run_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
run_dirs <- run_dirs[grepl("^run", basename(run_dirs))]

if (length(run_dirs) == 0) {
  stop("No run* folders found in: ", base_dir)
}

message("Found ", length(run_dirs), " run folders")

# Store per-run merged tables for later combine
all_run_tables <- list()

# =========================
# Process each run
# =========================
for (run_dir in run_dirs) {
  run_name <- basename(run_dir)
  run_tag  <- gsub("[^A-Za-z0-9]+", "_", run_name)

  message("\nProcessing: ", run_name)

  otu_file  <- file.path(run_dir, "08_otu_table.csv")
  taxa_file <- file.path(run_dir, "11_taxa.csv")

  if (!file.exists(otu_file)) {
    warning("Skipping ", run_name, ": missing 08_otu_table.csv")
    next
  }
  if (!file.exists(taxa_file)) {
    warning("Skipping ", run_name, ": missing 11_taxa.csv")
    next
  }

  # Read data
  otu_wide <- read_transpose_otu(otu_file)
  taxa_old <- fread(taxa_file, check.names = FALSE)

  # Make ASV values unique per run
  if ("ASV" %in% colnames(taxa_old)) {
    taxa_old[, ASV := paste0(ASV, "_", run_tag)]
  }

  tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  taxa_filled <- taxa_old %>%
    mutate(
      Kingdom = if_else(
        is.na(Kingdom) |
          str_trim(Kingdom) == "" |
          tolower(str_remove(Kingdom, "^[a-z]__")) %in% c("unassigned", "unclassified"),
        "d__Unassigned",
        Kingdom
      )
    ) %>%
    mutate(
      across(
        Phylum:Species,
        ~ {
          x <- replace_na(.x, "")
          core <- str_remove(x, "^[a-z]__")
          core <- str_trim(core)
          if_else(tolower(core) %in% c("", "unassigned", "unclassified"), "unclassified", core)
        }
      )
    ) %>%
    mutate(
      Kingdom = if_else(grepl("^[dpcofgs]__", Kingdom), Kingdom, paste0("d__", Kingdom)),
      Phylum  = if_else(grepl("^[dpcofgs]__", Phylum),  Phylum,  paste0("p__", Phylum)),
      Class   = if_else(grepl("^[dpcofgs]__", Class),   Class,   paste0("c__", Class)),
      Order   = if_else(grepl("^[dpcofgs]__", Order),   Order,   paste0("o__", Order)),
      Family  = if_else(grepl("^[dpcofgs]__", Family),  Family,  paste0("f__", Family)),
      Genus   = if_else(grepl("^[dpcofgs]__", Genus),   Genus,   paste0("g__", Genus)),
      Species = if_else(grepl("^[dpcofgs]__", Species), Species, paste0("s__", Species))
    )

  taxa_filled[, tax_cols] <- as.data.table(
    t(apply(taxa_filled[, tax_cols, with = FALSE], 1, fill_from_parent))
  )

  taxa <- taxa_filled

  if (!"Sequence" %in% colnames(taxa)) {
    stop("Sequence column not found in ", taxa_file)
  }

  # Merge using Sequence
  merged <- merge(
    taxa,
    otu_wide,
    by = "Sequence",
    all.y = TRUE,
    sort = FALSE
  )

  desired_meta <- c("ASV", "Sequence", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  existing_meta <- desired_meta[desired_meta %in% colnames(merged)]
  sample_cols <- setdiff(colnames(merged), existing_meta)

  merged <- merged[, c(existing_meta, sample_cols), with = FALSE]

  # Save per-run otu+taxa counts table
  out_counts <- file.path(run_dir, "12_otu_taxa_table.csv")
  fwrite(merged, out_counts, quote = TRUE)

  # Save per-run relative abundance table
  relabund <- make_relabund(
    merged,
    meta_cols = existing_meta
  )
  out_relabund <- file.path(run_dir, "13_otu_taxa_relabund.csv")
  fwrite(relabund, out_relabund, quote = TRUE)

  # Store for combined full join
  all_run_tables[[run_name]] <- copy(merged)
}

# =========================
# Combine all runs by Sequence (full join)
# =========================
if (length(all_run_tables) == 0) {
  stop("No run tables were successfully created.")
}

message("\nCombining all runs by Sequence...")

tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# ---- 1) Build a master taxonomy table by Sequence ----
taxonomy_tables <- lapply(all_run_tables, function(df) {
  unique(df[, c("Sequence", tax_cols), with = FALSE])
})

taxonomy_long <- rbindlist(taxonomy_tables, fill = TRUE, use.names = TRUE)

taxonomy_master <- taxonomy_long[
  ,
  lapply(.SD, function(x) {
    vals <- x[!is.na(x) & x != ""]
    if (length(vals) == 0) "" else vals[1]
  }),
  by = Sequence,
  .SDcols = tax_cols
]

# ---- 2) Build run-specific ASV + count tables ----
tables_for_merge <- lapply(names(all_run_tables), function(run_name) {
  df <- copy(all_run_tables[[run_name]])
  run_tag <- gsub("[^A-Za-z0-9]+", "_", run_name)

  asv_col <- paste0("ASV_", run_tag)
  df[, (asv_col) := ASV]
  df[, ASV := NULL]

  sample_cols <- setdiff(colnames(df), c("Sequence", tax_cols, asv_col))

  out <- unique(df[, c("Sequence", asv_col, sample_cols), with = FALSE])

  return(out)
})

# ---- 3) Full join all run-specific tables by Sequence ----
combined_counts_asv <- Reduce(
  function(x, y) merge(x, y, by = "Sequence", all = TRUE, sort = FALSE),
  tables_for_merge
)

# ---- 4) Join taxonomy once at the end ----
combined <- merge(
  taxonomy_master,
  combined_counts_asv,
  by = "Sequence",
  all = TRUE,
  sort = FALSE
)

# Fill missing per-run ASV columns with empty string
asv_cols <- grep("^ASV_", colnames(combined), value = TRUE)
for (col in asv_cols) {
  set(combined, which(is.na(combined[[col]])), col, "")
}

# Fill missing sample counts with 0
sample_cols <- setdiff(colnames(combined), c("Sequence", tax_cols, asv_cols))
for (s in sample_cols) {
  set(combined, which(is.na(combined[[s]])), s, 0)
}

# Add stable combined ID
combined[, combined_ASV := paste0("combined_ASV", seq_len(.N))]

# Reorder columns and drop Sequence from final combined outputs
combined <- combined[, c(
  "combined_ASV",
  tax_cols,
  asv_cols,
  sample_cols
), with = FALSE]

combined_counts_file <- file.path(base_dir, "12_combined_otu_taxa_table.csv")
fwrite(combined, combined_counts_file, quote = TRUE)

# Combined relative abundance
combined_relabund <- make_relabund(
  combined,
  meta_cols = c("combined_ASV", tax_cols, asv_cols)
)

combined_relabund_file <- file.path(base_dir, "13_combined_otu_taxa_relabund.csv")
fwrite(combined_relabund, combined_relabund_file, quote = TRUE)

message("\nDone.")
message("Per-run outputs written inside each run folder:")
message("  12_otu_taxa_table.csv")
message("  13_otu_taxa_relabund.csv")
message("\nCombined outputs written to:")
message("  ", combined_counts_file)
message("  ", combined_relabund_file)