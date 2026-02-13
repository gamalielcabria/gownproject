#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

stats_path <- "/home/gam/github/gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/01_SeqData/02_run200316/truncF270_truncR230.stats.tsv"
outdir <- "00_stats_hist"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- read stats (tab-delimited; skip comment lines like #q2:types) ----
st <- read.delim(
  stats_path,
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ---- normalize names ----
names(st) <- names(st) %>%
  stringr::str_trim() %>%
  stringr::str_replace_all("\\s+", "_") %>%
  stringr::str_replace_all("-", "_")

# Identify columns robustly
# Expect: sample_id + non_chimeric + percentage_of_input_non_chimeric
sample_col <- names(st)[1]

if (!("non_chimeric" %in% names(st))) {
  stop("Missing column 'non-chimeric' (normalized to non_chimeric). Found: ", paste(names(st), collapse=", "))
}
if (!("percentage_of_input_non_chimeric" %in% names(st))) {
  stop("Missing column 'percentage of input non-chimeric' (normalized to percentage_of_input_non_chimeric). Found: ", paste(names(st), collapse=", "))
}

# Check column count
if (ncol(st) < 9) {
  stop("Unexpected number of columns in stats file. Found ", ncol(st))
}

st2 <- tibble(
  sample_id = st[[1]],
  non_chimeric_reads = as.numeric(st[[8]]),
  pct_non_chimeric   = as.numeric(st[[9]])
) %>%
  filter(!is.na(non_chimeric_reads), !is.na(pct_non_chimeric))

base <- basename(stats_path) %>% stringr::str_replace("\\.stats\\.tsv$", "")

# Save extracted values
st2_ranked <- st2 %>%
  arrange(non_chimeric_reads) %>%
  mutate(rank_low_to_high = row_number())

write.table(
  st2_ranked,
  file = file.path(outdir, paste0(base, "_nonchim_reads_and_pct_ranked.tsv")),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---- Histogram: non-chimeric reads ----
p_reads <- ggplot(st2, aes(x = non_chimeric_reads)) +
  geom_histogram(bins = 30) +
  labs(
    title = paste0(base, ": non-chimeric reads per sample"),
    subtitle = paste0(
      "n=", nrow(st2),
      " | median=", scales::comma(round(median(st2$non_chimeric_reads))),
      " | mean=", scales::comma(round(mean(st2$non_chimeric_reads)))
    ),
    x = "Non-chimeric reads",
    y = "Sample count"
  ) +
  theme_bw()

ggsave(
  filename = file.path(outdir, paste0(base, "_hist_nonchim_reads.png")),
  plot = p_reads, width = 7, height = 5, dpi = 300
)

# ---- Histogram: % non-chimeric ----
p_pct <- ggplot(st2, aes(x = pct_non_chimeric)) +
  geom_histogram(bins = 30) +
  labs(
    title = paste0(base, ": % of input non-chimeric"),
    subtitle = paste0(
      "n=", nrow(st2),
      " | median=", sprintf("%.2f", median(st2$pct_non_chimeric)),
      " | mean=", sprintf("%.2f", mean(st2$pct_non_chimeric))
    ),
    x = "% of input non-chimeric",
    y = "Sample count"
  ) +
  theme_bw()

ggsave(
  filename = file.path(outdir, paste0(base, "_hist_pct_nonchim.png")),
  plot = p_pct, width = 7, height = 5, dpi = 300
)

# ============================================================
# APPENDED: log-scale read distribution + label extremes
# ============================================================

# pick extremes to label (min + max)
ext_reads <- st2 %>%
  slice(c(which.min(non_chimeric_reads), which.max(non_chimeric_reads))) %>%
  distinct(sample_id, .keep_all = TRUE)

ext_pct <- st2 %>%
  slice(c(which.min(pct_non_chimeric), which.max(pct_non_chimeric))) %>%
  distinct(sample_id, .keep_all = TRUE)

# ---- (A) Reads histogram on log10 x-axis + extreme labels ----
p_reads_log <- ggplot(st2, aes(x = non_chimeric_reads)) +
  geom_histogram(bins = 30) +
  scale_x_log10(labels = scales::comma) +
  geom_vline(xintercept = median(st2$non_chimeric_reads), linetype = 2) +
  geom_text(
    data = ext_reads,
    aes(x = non_chimeric_reads, y = 0, label = sample_id),
    inherit.aes = FALSE,
    angle = 90, vjust = -0.2, hjust = 0,
    size = 3
  ) +
  labs(
    title = paste0(base, ": non-chimeric reads per sample (log10 x-axis)"),
    subtitle = paste0(
      "n=", nrow(st2),
      " | median=", scales::comma(round(median(st2$non_chimeric_reads))),
      " | min=", scales::comma(round(min(st2$non_chimeric_reads))),
      " | max=", scales::comma(round(max(st2$non_chimeric_reads)))
    ),
    x = "Non-chimeric reads per sample (log10 scale)",
    y = "Sample count"
  ) +
  theme_bw()

ggsave(
  filename = file.path(outdir, paste0(base, "_hist_nonchim_reads_log10_labeled.png")),
  plot = p_reads_log, width = 7.5, height = 5.5, dpi = 300
)

# ---- (B) % non-chimeric histogram + extreme labels ----
p_pct_labeled <- ggplot(st2, aes(x = pct_non_chimeric)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = median(st2$pct_non_chimeric), linetype = 2) +
  geom_text(
    data = ext_pct,
    aes(x = pct_non_chimeric, y = 0, label = sample_id),
    inherit.aes = FALSE,
    angle = 90, vjust = -0.2, hjust = 0,
    size = 3
  ) +
  labs(
    title = paste0(base, ": % of input non-chimeric (labeled extremes)"),
    subtitle = paste0(
      "n=", nrow(st2),
      " | median=", sprintf("%.2f", median(st2$pct_non_chimeric)),
      " | min=", sprintf("%.2f", min(st2$pct_non_chimeric)),
      " | max=", sprintf("%.2f", max(st2$pct_non_chimeric))
    ),
    x = "% of input non-chimeric",
    y = "Sample count"
  ) +
  theme_bw()

ggsave(
  filename = file.path(outdir, paste0(base, "_hist_pct_nonchim_labeled.png")),
  plot = p_pct_labeled, width = 7.5, height = 5.5, dpi = 300
)

message("Done. Outputs in: ", normalizePath(outdir))
