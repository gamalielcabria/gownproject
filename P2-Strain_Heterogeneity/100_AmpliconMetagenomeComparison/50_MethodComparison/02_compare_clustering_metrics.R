library(tidyverse)
library(aricode)
library(vegan)
library(dynamicTreeCut)
library(gghalves)
library(ggsci)

# Import pre-computed microeco objects
rdata_path <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/50_MethodComparison/meco_objects.rds"
meco_objects <- readRDS(rdata_path)
list2env(meco_objects, envir = .GlobalEnv)

figure_path <- "/home/gam/github/gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/00_Figures"

#########################################
# COMPARE TAXA-LEVEL CLUSTERING METRICS #
#########################################
# This section computes clustering metrics (ARI, NMI) at the genus level
# for each method's OTU_table compared to the amplicon-based OTU_table (VSEARCH).

## Function to compute comparative metrics by Taxa
collapse_to_taxa <- function(otu, tax, taxa_col = "Genus") {
  stopifnot(all(rownames(otu) %in% rownames(tax)))

  g <- tax[rownames(otu), taxa_col]
  g <- as.character(g)

  # clean genus labels
  g <- sub("^[a-z]__", "", g)
  g[is.na(g) | g == "" | tolower(g) == "unclassified"] <- "Unclassified"

  # sum features to genus
  # result: genus x samples numeric matrix
  genus_mat <- rowsum(as.matrix(otu), group = g, reorder = TRUE)

  # drop genus with all zeros (optional)
  genus_mat <- genus_mat[rowSums(genus_mat) > 0, , drop = FALSE]
  genus_mat
}

taxa_richness_per_sample <- function(taxa_mat, min_abund = 1) {
  # taxa_mat: taxa x samples
  colSums(taxa_mat >= min_abund)
}

## Prepare OTU and Tax inputs for comparison 
### Choose rank
#taxon_level <- "Genus"   
min_abund   <- 0.001         # presence threshold after collapsing (counts); keep 1 unless you have fractional abund

### Build long dataframe: sample, method, richness
taxon_level <- "Genus" # Options are: c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rich_long_genus <- imap_dfr(meco_objects, function(meco, method_name) {
  otu <- meco$otu_table
  tax <- meco$tax_table

  collapsed <- collapse_to_taxa(otu, tax, taxa_col = taxon_level)
  rich <- taxa_richness_per_sample(collapsed, min_abund = min_abund)

  tibble(
    sample = names(rich),
    method = method_name,
    richness = as.numeric(rich),
    rank = taxon_level
  )
})

taxon_level <- "Order" # Options are: c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rich_long_order <- imap_dfr(meco_objects, function(meco, method_name) {
  otu <- meco$otu_table
  tax <- meco$tax_table

  collapsed <- collapse_to_taxa(otu, tax, taxa_col = taxon_level)
  rich <- taxa_richness_per_sample(collapsed, min_abund = min_abund)

  tibble(
    sample = names(rich),
    method = method_name,
    richness = as.numeric(rich),
    rank = taxon_level
  )
})

taxon_level <- "Phylum" # Options are: c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rich_long_phylum <- imap_dfr(meco_objects, function(meco, method_name) {
  otu <- meco$otu_table
  tax <- meco$tax_table

  collapsed <- collapse_to_taxa(otu, tax, taxa_col = taxon_level)
  rich <- taxa_richness_per_sample(collapsed, min_abund = min_abund)

  tibble(
    sample = names(rich),
    method = method_name,
    richness = as.numeric(rich),
    rank = taxon_level
  )
})

rich_long <- bind_rows(rich_long_genus, rich_long_order, rich_long_phylum)

### Keeping only samples shared by *all* methods (recommended for fair comparison)
shared_samples <- rich_long %>%
  count(sample) %>%
  filter(n == ( n_distinct(rich_long$method) * n_distinct(rich_long$rank))) %>%
  pull(sample)

rich_long_shared <- rich_long %>%
  filter(sample %in% shared_samples) %>%
  separate(method, into = c("Method", "OTU", "Database"), sep = "_", extra = "drop", remove = FALSE) %>%
  mutate(Methods = paste0(str_to_sentence(Method), "(", str_to_upper(Database),")")) %>%
  select(-OTU, -Method) %>%
  mutate(Sequencing = ifelse(grepl("dada2|vsearch", Methods, ignore.case = TRUE), "Amplicon", "Methods"))

## Plot richness comparison
### Prepare method order
method_order <- rich_long_shared %>%
  distinct(Sequencing, Methods) %>%                # one row per method
  mutate( Sequencing = factor( Sequencing, levels = c("Amplicon", "Metagenome") )  ) %>%
  arrange(Sequencing, Methods) %>%                  # 1) sequencing, 2) alpha
  pull(Methods)

### Barplot
taxon_level <- paste(sort(unique(na.omit(rich_long_shared$rank))), collapse = ", ")

plot_richness_boxplot <- ggplot(
    rich_long_shared %>% mutate(Methods = factor(Methods, levels = method_order)), 
    aes(x = Methods, y = richness, fill = Database)
    ) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.1) +
  geom_violin(width = 1, alpha = 0.4) +
  geom_line(aes(group=sample),alpha = 0.2) +
  geom_point(alpha = 0.7, size = 1.2) +
  #geom_jitter(width = 0.15, height = 0, alpha = 0.6) +
  scale_fill_jco(
    name = "Database",
    labels = c("GTDB", "SGB", "SILVA")    
  ) +
  labs(
    title = paste("Per-sample taxa richness"), #at", taxon_level),
    x = "Method",
    y = "Taxa Richness"
  ) +
  facet_wrap(~rank, scales = "free_y", nrow = 3) +
  theme_minimal() +
  theme(
    #aspect.ratio = 3/2,
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(10, 10, 10, 10),
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top",
    legend.direction = "horizontal"
    )

# ggsave(
#   file = file.path(figure_path, "/richness_boxplot.png"),
#   plot_richness_boxplot,
#   width = 5,
#   height = 12,
#   dpi = 300
# )

# saveRDS(
#   plot_richness_boxplot,
#   file = file.path(figure_path, "/Figure_SX_taxa_richness_comparison_boxplot.rds")
# )


# Comparing Alpha Diversity between Methods

# alpha metrics from an OTU table (taxa x samples OR samples x taxa)
alpha_measures <- c("Shannon", "Simpson", "Pielou", "Coverage")#, "Chao1", "ACE", "Coverage")
# Add Fisher if it works for your data; sometimes it can fail depending on counts.
# alpha_measures <- c(alpha_measures, "Fisher")

### SANITY CHECK: any negative values in OTU tables?
neg_check <- imap_dfr(meco_objects, function(meco, method_name) {
  x <- as.matrix(meco$otu_table)
  tibble(
    method = method_name,
    min_value = suppressWarnings(min(x, na.rm = TRUE)),
    n_negative = sum(x < 0, na.rm = TRUE),
    n_na = sum(is.na(x))
  )
})
sanity_check <- neg_check %>% arrange(desc(n_negative), min_value)
sanity_check

### Compute alpha diversity for each method
alpha_long <- imap_dfr(meco_objects, function(meco, method_name) {

  # --- sanitize otu_table for alpha diversity ---
  otu <- as.matrix(meco$otu_table)
  otu[is.na(otu)] <- 0
  otu[otu < 0] <- 0
  meco$otu_table <- otu

  # compute alpha via microeco
  meco$cal_alphadiv(measures = alpha_measures, PD = FALSE)

  meco$alpha_diversity %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
    mutate(method = method_name)
})

### Keep only shared samples
shared_samples_alpha <- alpha_long %>%
  count(sample, metric) %>%
  group_by(metric) %>%
  filter(n == n_distinct(alpha_long$method)) %>%
  ungroup() %>%
  distinct(sample) %>%
  pull(sample)

alpha_long_shared <- alpha_long %>%
  filter(sample %in% shared_samples_alpha) %>%
  separate(method, into = c("Method", "OTU", "Database"), sep = "_", extra = "drop", remove = FALSE) %>%
  mutate(
    Methods = paste0(str_to_sentence(Method), "(", str_to_upper(Database), ")"),
    Sequencing = ifelse(grepl("dada2|vsearch", method, ignore.case = TRUE), "Amplicon", "Metagenome")
  ) %>%
  select(-OTU, -Method)

# Prep statistical comparisons
### Setup: !!!!!!!!!!!! EDIT THIS !!!!!!!!!!!!!!!!
ref_method <- "Vsearch(GTDB)"

compare_methods <- c(
  "Singlem(GTDB)",
  "Mp4(GTDB)",
  "Dada2(GTDB)"
)

### Build paired dataframe for statistical tests
alpha_paired <- alpha_long_shared %>%
  filter(Methods %in% c(ref_method, compare_methods)) %>%
  select(sample, metric, Methods, value) %>%
  pivot_wider(
    names_from = Methods,
    values_from = value
  )

### Perform paired statistical tests
stat_results <- map_dfr(compare_methods, function(m) {
  alpha_paired %>%
    group_by(metric) %>%
    summarise(
      method = m,
      p_value = wilcox.test(
        .data[[m]],
        .data[[ref_method]],
        paired = TRUE,
        exact = FALSE
      )$p.value,
      .groups = "drop"
    )
})

stat_results <- stat_results %>%
  group_by(metric) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

stat_annot <- alpha_long_shared %>%
  filter(Methods %in% compare_methods) %>%
  group_by(metric, Methods) %>%
  summarise(y_pos = max(value, na.rm = TRUE) * 1.08, .groups = "drop") %>%
  left_join(
    stat_results,
    by = c("metric", "Methods" = "method")
  ) %>%
  mutate(
    label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

# Plot alpha diversity comparison
### Prepare method order
method_order_alpha <- alpha_long_shared %>%
  distinct(Sequencing, Methods) %>%
  mutate(Sequencing = factor(Sequencing, levels = c("Amplicon", "Metagenome"))) %>%
  arrange(Sequencing, Methods) %>%
  pull(Methods)

### Main Plot: Alpha Diversity Comparison
plot_alpha <- ggplot(
  alpha_long_shared %>% mutate(Methods = factor(Methods, levels = method_order_alpha)),
  aes(x = Methods, y = value, fill = Database)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  #geom_violin(alpha = 0.35, width = 1) +
  #geom_line(aes(group = sample), alpha = 0.18) +
  geom_point(alpha = 0.7, size = 1.1) +
  scale_fill_jco(
    name = "Database",
    labels = c("GTDB", "SGB", "SILVA")    
  ) +
  facet_wrap(~ metric, scales = "free_y", nrow = 2) +
  labs(
    title = "Per-sample alpha diversity comparison across methods",
    x = "Method",
    y = "Alpha diversity"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1, "lines"),
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top",
    legend.direction = "horizontal"
  )

plot_alpha

### Add Statistical Annotations
plot_alpha_stat <- plot_alpha +
geom_text(
    data = stat_annot,
    aes(
      x = Methods,
      y = y_pos,
      label = label
    ),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  )

# ggsave(
#   file = file.path(figure_path, "Figure_SX_alpha_diversity_comparison_boxplot.png"),
#   plot_alpha_stat,
#   width = 6,
#   height = 8,
#   dpi = 300
# )

# saveRDS(
#   plot_alpha_stat,
#   file = file.path(figure_path, "Figure_SX_alpha_diversity_comparison_boxplot.rds")
# )


###############################################################
# Beta Diversity Comparison can be added similarly if needed. #
###############################################################

## Function to compute and extract Bray–Curtis dist for one microeco object
get_bray_dist <- function(meco, samples_keep = NULL) {
  # Ensure non-negative (Bray requires non-negative abundances)
  otu <- as.matrix(meco$otu_table)
  otu[is.na(otu)] <- 0
  if (any(otu < 0)) otu[otu < 0] <- 0
  meco$otu_table <- otu

  # Compute Bray–Curtis inside microeco
  meco$cal_betadiv(method = "bray")  # stores into meco$beta_diversity :contentReference[oaicite:2]{index=2}

  # microeco stores distance matrices in a list; Bray is typically under $bray
  bray <- meco$beta_diversity$bray

  # keep shared samples + enforce same order
  if (!is.null(samples_keep)) {
    smp <- intersect(samples_keep, rownames(as.matrix(bray)))
    bray <- as.matrix(bray)[smp, smp, drop = FALSE]
    bray <- as.dist(bray)
  }
  bray
}

## Shared samples across all methods for beta diversity
all_samples_by_method <- imap(meco_objects, ~ colnames(.x$otu_table))
shared_samples <- Reduce(intersect, all_samples_by_method)

## Compute distances for each method
bray_by_method <- imap(meco_objects, ~ get_bray_dist(.x, samples_keep = shared_samples))

# Compare two methods pairwise via Mantel test
compare_methods_mantel <- function(bray_list, methodA, methodB,
                                   cor_method = "spearman", permutations = 9999) {
  dA <- bray_list[[methodA]]
  dB <- bray_list[[methodB]]

  # enforce identical labels/order
  sA <- attr(dA, "Labels")
  sB <- attr(dB, "Labels")
  s  <- intersect(sA, sB)

  dA <- as.dist(as.matrix(dA)[s, s])
  dB <- as.dist(as.matrix(dB)[s, s])

  vegan::mantel(dA, dB, method = cor_method, permutations = permutations)
}

### EDIT THIS!!!
mantel_res <- compare_methods_mantel(bray_by_method, "vsearch_OTU_gtdb_meco", "singlem_OTU_gtdb_meco")
mantel_res

### PLOT Mantel results across multiple method comparisons
dist_scatter_df <- function(bray_list, methodA, methodB) {
  dA <- bray_list[[methodA]]
  dB <- bray_list[[methodB]]

  s  <- intersect(attr(dA, "Labels"), attr(dB, "Labels"))
  mA <- as.matrix(dA)[s, s]
  mB <- as.matrix(dB)[s, s]

  idx <- upper.tri(mA, diag = FALSE)

  tibble(
    pair = paste(rep(s, each = length(s))[idx], rep(s, times = length(s))[idx], sep = " vs "),
    dist_A = mA[idx],
    dist_B = mB[idx]
  )
}
method_a <- "vsearch_OTU_gtdb_meco"
method_b <- "singlem_OTU_gtdb_meco"
df_scatter <- dist_scatter_df(bray_by_method, method_a, method_b)

pretty_method <- function(x) {
  tibble(raw = x) %>%
    separate(raw, into = c("method", "OTU", "db", "rest"), sep = "_", fill = "right") %>%
    transmute(label = paste0(toupper(method), "(", toupper(db), ")")) %>%
    pull(label)
}

cor_res <- cor.test(
  df_scatter$dist_A,
  df_scatter$dist_B,
  method = "spearman",
  use = "complete.obs"
)

rho <- unname(cor_res$estimate)
pval <- cor_res$p.value

lims <- range(c(df_scatter$dist_A, df_scatter$dist_B), na.rm = TRUE)

plot_beta_bray_scatter <- ggplot(df_scatter, aes(dist_A, dist_B)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  coord_equal(xlim = lims, ylim = lims) +
  annotate(
    "text",
    x = lims[1] + 0.05 * diff(lims),
    y = lims[2] - 0.05 * diff(lims),
    hjust = 0,
    vjust = 1,
    label = sprintf("Spearman \u03C1 = %.2f\np = %.2g", rho, pval),
    size = 4,
    fontface = "bold"
  ) +
  labs(
    x = paste0("Bray–Curtis ", pretty_method(method_a)),
    y = paste0("Bray–Curtis ", pretty_method(method_b)),
    title = paste0("Pairwise sample distances:\n",pretty_method(method_a)," vs ",pretty_method(method_b))
  ) +
  theme_minimal()

plot_beta_bray_scatter

# ggsave(
#   file = file.path(figure_path, "Figure_SX_beta_diversity_bray_scatter.png"),
#   plot_beta_bray_scatter,
#   width = 6,
#   height = 5,
#   dpi = 300
# )

# saveRDS(
#   plot_beta_bray_scatter,
#   file = file.path(figure_path, "Figure_SX_beta_diversity_bray_scatter.rds")
# )

# Comparison of multiple samples to a reference:
ref <- "vsearch_OTU_gtdb_meco"
others <- setdiff(names(bray_by_method), ref)

mantel_table <- map_dfr(others, function(m) {
  res <- compare_methods_mantel(bray_by_method, ref, m, cor_method = "spearman", permutations = 9999)
  tibble(
    ref = ref,
    method = m,
    mantel_r = unname(res$statistic),
    p_value = res$signif
  )
}) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

mantel_table

# write.csv(   mantel_table,   row.names = FALSE
#   file = file.path(figure_path, "Table_SX_beta_diversity_mantel_results.csv"),
# )
