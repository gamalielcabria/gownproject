#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})
setwd("/home/gam/github/gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/01_SeqData/02_run200316")

# ----------------------------
# Inputs
# ----------------------------
infile <- "trunc_sweep_summary.tsv"   # your summary TSV
outdir <- "00_trunc_sweep_viz"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Load + parse truncF/truncR
# ----------------------------
df <- readr::read_tsv(infile, show_col_types = FALSE) %>%
  separate(truncF_truncR, into = c("truncF", "truncR"), sep = "_", convert = TRUE) %>%
  mutate(
    across(c(mean_passed_filter, mean_merged, mean_non_chimeric, repseq_count), as.numeric)
  )

# ----------------------------
# Robust scoring
# - We want "best combination" based on metrics.
# - Use ranks (robust to scale differences).
# - Higher is better for all metrics here.
# - Overall score = average of rank-percentiles.
# ----------------------------
rank01 <- function(x) {
  # rank percentile in [0,1], higher is better
  r <- rank(x, ties.method = "average", na.last = "keep")
  (r - 1) / (sum(!is.na(x)) - 1)
}

df_scored <- df %>%
  mutate(
    r_pass = rank01(mean_passed_filter),
    r_merge = rank01(mean_merged),
    r_nonch = rank01(mean_non_chimeric),
    r_rep   = rank01(log10(repseq_count)),
    score   = (r_pass + r_merge + r_nonch + r_rep) / 4,
    combo   = paste0(truncF, "_", truncR)
  ) %>%
  arrange(desc(score))

best <- df_scored %>% slice(1)

# Save ranked table
readr::write_tsv(df_scored %>% select(combo, truncF, truncR, mean_passed_filter, mean_merged, mean_non_chimeric, repseq_count, score),
                 file.path(outdir, "ranked_combinations.tsv"))

# ----------------------------
# 1) Heatmaps of raw metrics with "best" outlined
# ----------------------------
df_long <- df_scored %>%
  select(truncF, truncR, combo, mean_passed_filter, mean_merged, mean_non_chimeric, repseq_count) %>%
  pivot_longer(
    cols = c(mean_passed_filter, mean_merged, mean_non_chimeric, repseq_count),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
      mean_passed_filter = "Mean % passed filter",
      mean_merged = "Mean % merged",
      mean_non_chimeric = "Mean % non-chimeric",
      repseq_count = "Rep-seq count"
    )
  )

# data frame for best-tile outline
best_tile <- best %>% transmute(truncF, truncR)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggnewscale)   # install.packages("ggnewscale") if needed
})

# df_long: columns truncF, truncR, metric, value
# best_tile: columns truncF, truncR

pct_metrics <- c("Mean % passed filter", "Mean % merged", "Mean % non-chimeric")
rep_metric  <- "Rep-seq count"

df_long2 <- df_long %>%
  mutate(
    metric = factor(metric, levels = c(pct_metrics, rep_metric)),
    label = case_when(
      metric %in% pct_metrics ~ sprintf("%.1f", value),
      metric == rep_metric    ~ scales::comma(round(value)),
      TRUE ~ as.character(value)
    ),
    value_log10 = if_else(metric == rep_metric, log10(value), NA_real_)
  )

p_heat <- ggplot() +
  # ---- percentage facets (normal fill scale) ----
  geom_tile(
    data = df_long2 %>% filter(metric %in% pct_metrics),
    aes(x = truncR, y = truncF, fill = value)
  ) +
  geom_text(
    data = df_long2 %>% filter(metric %in% pct_metrics),
    aes(x = truncR, y = truncF, label = label),
    size = 3
  ) +
  ggnewscale::new_scale_fill() +
  # ---- repseq facet (log10 fill scale) ----
  geom_tile(
    data = df_long2 %>% filter(metric == rep_metric),
    aes(x = truncR, y = truncF, fill = value_log10)
  ) +
  geom_text(
    data = df_long2 %>% filter(metric == rep_metric),
    aes(x = truncR, y = truncF, label = label),
    size = 3
  ) +
  # outline best combo on all facets
  geom_point(
    data = best_tile,
    aes(x = truncR, y = truncF),
    inherit.aes = FALSE,
    size = 6,
    stroke = 1.5,
    shape = 21,
    fill = NA
  ) +
  facet_wrap(~ metric, scales = "free") +
  labs(
    title = "Truncation sweep metrics (values in tiles)",
    subtitle = paste0("Best by composite score: ", best$combo,
                      " | score=", sprintf("%.3f", best$score)),
    x = "truncR",
    y = "truncF"
  ) +
  # keep legends clear (no explicit colors set)
  guides(
    fill = guide_colorbar(title = "Value (log10 for repseq facet)")
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_heat

ggsave(file.path(outdir, "heatmaps_metrics_best_outlined.png"), p_heat, width = 12, height = 8, dpi = 300)

# ----------------------------
# 2) Heatmap of overall score (easy “best combo” view)
# ----------------------------
p_score <- ggplot(df_scored, aes(x = truncR, y = truncF, fill = score)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", score)), size = 3) +
  geom_point(
    data = best_tile,
    aes(x = truncR, y = truncF),
    inherit.aes = FALSE,
    size = 6,
    stroke = 1.5,
    shape = 21,
    fill = NA
  ) +
  labs(
    title = "Composite score heatmap (rank-based; higher = better)",
    subtitle = paste0("Best: ", best$combo, " (score=", sprintf("%.3f", best$score), ")"),
    x = "truncR",
    y = "truncF",
    fill = "Score"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave(file.path(outdir, "heatmap_composite_score.png"), p_score, width = 7.5, height = 6, dpi = 300)

# ----------------------------
# 3) Top-N bar plot of combos by score
# ----------------------------
top_n <- 10

p_top <- df_scored %>%
  slice_head(n = top_n) %>%
  mutate(combo = fct_reorder(combo, score)) %>%
  ggplot(aes(x = combo, y = score)) +
  geom_col() +
  coord_flip() +
  labs(
    title = paste0("Top ", top_n, " truncF_truncR combinations by composite score"),
    x = "Combination (truncF_truncR)",
    y = "Composite score"
  ) +
  theme_bw()

ggsave(file.path(outdir, "top_combos_by_score.png"), p_top, width = 7.5, height = 5.5, dpi = 300)

# ----------------------------
# 4) Pareto front: maximize non-chimeric % and repseq_count
# (shows best trade-offs; no weighting)
# ----------------------------
pareto <- df_scored %>%
  arrange(desc(mean_non_chimeric), desc(repseq_count)) %>%
  mutate(
    max_rep_so_far = cummax(repseq_count),
    on_pareto = repseq_count >= max_rep_so_far
  )

p_pareto <- ggplot(pareto, aes(x = mean_non_chimeric, y = repseq_count)) +
  geom_point(aes(shape = on_pareto), size = 2) +
  labs(
    title = "Pareto view: maximize mean % non-chimeric and rep-seq count",
    subtitle = paste0("Best-by-score outlined elsewhere: ", best$combo),
    x = "Mean % non-chimeric",
    y = "Rep-seq count",
    shape = "On Pareto front"
  ) +
  theme_bw()

ggsave(file.path(outdir, "pareto_nonchim_vs_repseq.png"), p_pareto, width = 7.5, height = 5.5, dpi = 300)

# ----------------------------
# Print best combo to stdout
# ----------------------------
cat(
  "Best combination by composite rank score:\n",
  best$combo, "\n",
  "truncF=", best$truncF, " truncR=", best$truncR, "\n",
  "mean_passed_filter=", best$mean_passed_filter, "\n",
  "mean_merged=", best$mean_merged, "\n",
  "mean_non_chimeric=", best$mean_non_chimeric, "\n",
  "repseq_count=", best$repseq_count, "\n",
  "score=", best$score, "\n",
  sep = ""
)

message("Done. Outputs in: ", outdir)
