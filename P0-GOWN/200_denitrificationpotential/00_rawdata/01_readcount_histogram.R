#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# -----------------------------
# Inputs / outputs
# -----------------------------
infile  <- "gownproject/P0-GOWN/200_denitrificationpotential/00_rawdata/00_readcount_2024-2025.csv"
output_tag <- basename(infile) %>% str_remove("\\.csv$")
outdir <- "gownproject/P0-GOWN/200_denitrificationpotential/00_rawdata/"
out_tsv <- paste0(outdir, output_tag, "_with_metadata.tsv")
out_png <- paste0(outdir, output_tag, "_read_histogram.png")
out_csv_low <- paste0(outdir, output_tag, "_lowest10_samples.csv")
out_csv_high <- paste0(outdir, output_tag, "_highest10_samples.csv")

reads <- read_csv(infile, show_col_types = FALSE)

# -----------------------------
# Parse metadata from filename
# -----------------------------
reads_parsed <- reads %>%
  mutate(
    # Extract 2-digit year before "GOWN" (e.g., "25" in "GLC_25GOWN...")
    Year_short = str_extract(R1_file, "(?<=GLC_)[0-9]{2}(?=GOWN)"),
    Year       = as.integer(paste0("20", Year_short)),

    # Project token (e.g., "GOWN" in "GLC_25GOWN...")
    Project    = str_extract(R1_file, "(?<=GLC_[0-9]{2})[A-Z]+"),

    # Well (captures GW265, GW265NG1, GW3010, etc.)
    Well       = str_extract(R1_file, "GW[0-9A-Z]+"),

    # Genome code (captures G01, P56, L17, etc. between hyphens)
    GenomeCode = str_extract(R1_file, "(?<=-)[A-Z][0-9]{2}(?=-)"),

    # Convenience: reads in millions
    Read_M     = R1_reads / 1e6
  )

# -----------------------------
# Basic QC checks
# -----------------------------
# If R1 and R2 differ, flag them
reads_parsed <- reads_parsed %>%
  mutate(
    ReadDiff = R1_reads - R2_reads,
    PairOK   = if_else(ReadDiff == 0, TRUE, FALSE)
  )

message("Rows: ", nrow(reads_parsed))
message("Pairs with R1!=R2: ", sum(!reads_parsed$PairOK, na.rm = TRUE))
message("Missing parsed Year: ", sum(is.na(reads_parsed$Year)))
message("Missing parsed Project: ", sum(is.na(reads_parsed$Project)))
message("Missing parsed Well: ", sum(is.na(reads_parsed$Well)))
message("Missing parsed GenomeCode: ", sum(is.na(reads_parsed$GenomeCode)))

# -----------------------------
# Save table with metadata
# -----------------------------
# reads_parsed %>%
#   select(R1_file, R2_file, R1_reads, R2_reads, Match, Year, Project, Well, GenomeCode, Read_M, PairOK) %>%
#   write_tsv(out_tsv)

# message("Wrote metadata table: ", out_tsv)

# -----------------------------
# Plot: histogram of reads per sample
# -----------------------------
p <- ggplot(reads_parsed, aes(x = Read_M)) +
  geom_histogram(bins = 25, color = "black", fill = "#0073C2", linewidth = 0.35) +
  geom_vline(xintercept = median(reads_parsed$Read_M, na.rm = TRUE),
             linetype = "dashed", linewidth = 0.9) +
  labs(
    title = "Read Count Distribution per Sample (R1)",
    subtitle = paste0("n = ", nrow(reads_parsed),
                      " | median = ", round(median(reads_parsed$Read_M, na.rm = TRUE), 2), "M reads"),
    x = "Reads (Millions)",
    y = "Number of Samples"
  ) +
  theme_classic(base_size = 14)
p
ggsave(out_png, p, width = 8, height = 6, dpi = 300)
message("Wrote histogram: ", out_png)

# -----------------------------
# Extra: show lowest 10 samples (console)
# -----------------------------
message("\nLowest 10 samples by R1_reads:")
lowest10 <- reads_parsed %>%
  arrange(R1_reads) %>%
  transmute(
    Year, Project, Well, GenomeCode,
    R1_reads = comma(R1_reads),
    R1_M = round(Read_M, 2),
    R1_file
  ) %>%
  head(20) %>%
  print(n = 20)
lowest10

highest10 <- reads_parsed %>%
  arrange(desc(R1_reads)) %>%
  transmute(            
    Year, Project, Well, GenomeCode,
    R1_reads = comma(R1_reads),
    R1_M = round(Read_M, 2),
    R1_file
  ) %>%
  head(10) %>%
  print(n = 10)
highest10

write.csv(lowest10, out_csv_low, row.names = FALSE)
write.csv(highest10, out_csv_high, row.names = FALSE)