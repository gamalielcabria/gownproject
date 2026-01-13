# Create diversity
library(tidyverse)
library(microeco)

# Paths of combined otu tax files

vsearch_OTU_silva <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/20_qiime2_otu/silva_bact_arch_combined_qiimeOTU.csv"
vsearch_OTU_gtdb <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/20_qiime2_otu/gtdb_bact_arch_combined_qiimeOTU.csv"
dada2_ASV_gtdb <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv/gtdb_bact_arch_combined_qiimeASV.csv"
dada2_ASV_silva <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv/silva_bact_arch_combined_qiimeASV.csv"
mp4_OTU_sgb <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/40_metaphlan4_otu/mp4-species_bact_arch_combined_otu.csv"
mp4_OTU_gtdb <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/40_metaphlan4_otu/mp4-species_bact_arch_combined_relabgtdb.csv"
singlem_OTU_gtdb <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/11_SingleM_genes/gtdb_bact_arch_combined_singlem-sum-OTU.csv"

metadata_path <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/10_GenusAggregate_Sum/sample_table_Meta.csv"

# Read files and separate into otu and tax tables
## Read and Clean metadata
metadata_df <- read.csv(metadata_path, stringsAsFactors = FALSE) %>%
    separate(X, sep = "-", into = c("Sample_Code", "Treatment", "GOWN_Well"), extra = "drop") %>%
    select(Metagenome_Code, GOWN_Well, Year, Sample_Code) %>% 
    filter(!(Metagenome_Code %in% c("E14","E19","E21","E23") ) ) %>%
    column_to_rownames("Metagenome_Code")

meta_samples <- rownames(metadata_df)

## Read and Break into component taxa

tax_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")


### VSEARCH SILVA
vsearch_OTU_silva_df <- read.csv(vsearch_OTU_silva, stringsAsFactors = FALSE) %>%
    select(-X) %>%
    column_to_rownames("OTUs") %>%
    select(any_of(union(meta_samples,tax_cols)))

sample_cols <- setdiff(colnames(vsearch_OTU_silva_df), tax_cols)

vsearch_OTU_filtered <- vsearch_OTU_silva_df %>%
  rowwise() %>%
  mutate(total_count = sum(c_across(all_of(sample_cols)))) %>%
  ungroup() %>%
  filter(total_count > 2) %>%     # removes 0, 1, and 2
  select(-total_count)

vsearch_silva_OTU <- vsearch_OTU_filtered %>%
    select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
    as.data.frame()

vsearch_silva_TAX <- vsearch_OTU_filtered %>%
    select(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    as.data.frame()

vsearch_OTU_silva_meco <- microtable$new( otu_table = vsearch_silva_OTU, tax_table = vsearch_silva_TAX, sample_table = metadata_df)


### VSEARCH GTDB
vsearch_OTU_gtdb_df <- read.csv(vsearch_OTU_gtdb, stringsAsFactors = FALSE) %>%
    select(-X) %>%
    column_to_rownames("OTUs") %>%
    select(any_of(union(meta_samples,tax_cols)))

sample_cols <- setdiff(colnames(vsearch_OTU_gtdb_df), tax_cols)

vsearch_OTU_filtered <- vsearch_OTU_gtdb_df %>%
  rowwise() %>%
  mutate(total_count = sum(c_across(all_of(sample_cols)))) %>%
  ungroup() %>%
  filter(total_count > 2) %>%     # removes 0, 1, and 2
  select(-total_count)

vsearch_gtdb_OTU <- vsearch_OTU_filtered %>%
    select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
    as.data.frame()

vsearch_gtdb_TAX <- vsearch_OTU_filtered %>%
    select(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    as.data.frame()

vsearch_OTU_gtdb_meco <- microtable$new( otu_table = vsearch_gtdb_OTU, tax_table = vsearch_gtdb_TAX, sample_table = metadata_df)

### Dada2 ASV GTDB
dada2_ASV_gtdb_df <- read.csv(dada2_ASV_gtdb, stringsAsFactors = FALSE) %>%
    select(-X) %>%
    column_to_rownames("OTUs") %>%
    select(any_of(union(meta_samples,tax_cols)))

sample_cols <- setdiff(colnames(dada2_ASV_gtdb_df), tax_cols)

dada2_ASV_gtdb_filtered <- dada2_ASV_gtdb_df %>%
  rowwise() %>%
  mutate(total_count = sum(c_across(all_of(sample_cols)))) %>%
  ungroup() %>%
  filter(total_count > 2) %>%     # removes 0, 1, and 2
  select(-total_count)

dada2_ASV_gtdb_OTU <- dada2_ASV_gtdb_filtered %>%
    select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
    as.data.frame()

dada2_ASV_gtdb_TAX <- dada2_ASV_gtdb_filtered %>%
    select(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    as.data.frame()

dada2_ASV_gtdb_meco <- microtable$new( otu_table = dada2_ASV_gtdb_OTU, tax_table = dada2_ASV_gtdb_TAX, sample_table = metadata_df)

### DADA2 ASV Silva
dada2_ASV_silva_df <- read.csv(dada2_ASV_silva, stringsAsFactors = FALSE) %>%
    select(-X) %>%
    column_to_rownames("OTUs") %>%
    select(any_of(union(meta_samples,tax_cols)))

sample_cols <- setdiff(colnames(dada2_ASV_silva_df), tax_cols)

dada2_ASV_silva_filtered <- dada2_ASV_silva_df %>%
  rowwise() %>%
  mutate(total_count = sum(c_across(all_of(sample_cols)))) %>%
  ungroup() %>%
  filter(total_count > 2) %>%     # removes 0, 1, and 2
  select(-total_count)

dada2_ASV_silva_OTU <- dada2_ASV_silva_filtered %>%
    select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
    as.data.frame()

dada2_ASV_silva_TAX <- dada2_ASV_silva_filtered %>%
    select(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    as.data.frame()

dada2_ASV_silva_meco <- microtable$new( otu_table = dada2_ASV_silva_OTU, tax_table = dada2_ASV_silva_TAX, sample_table = metadata_df)


### MP4 OTU SGB
mp4_OTU_sgb_df <- read.csv(mp4_OTU_sgb, stringsAsFactors = FALSE) %>%
    select(-X) %>%
    column_to_rownames("OTUs") %>%
    select(any_of(union(meta_samples,tax_cols)))

sample_cols <- setdiff(colnames(mp4_OTU_sgb_df), tax_cols)

mp4_OTU_sgb_filtered <- mp4_OTU_sgb_df %>%
  rowwise() %>%
  mutate(total_count = sum(c_across(all_of(sample_cols)))) %>%
  ungroup() %>%
  filter(total_count > 2) %>%     # removes 0, 1, and 2
  select(-total_count)

mp4_OTU_sgb_OTU <- mp4_OTU_sgb_filtered %>%
    select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
    as.data.frame()

mp4_OTU_sgb_TAX <- mp4_OTU_sgb_filtered %>%
    select(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    as.data.frame()

mp4_OTU_sgb_meco <- microtable$new( otu_table = mp4_OTU_sgb_OTU, tax_table = mp4_OTU_sgb_TAX, sample_table = metadata_df)

### MP4 OTU GTDB
mp4_OTU_gtdb_df <- read.csv(mp4_OTU_gtdb, stringsAsFactors = FALSE) %>%
    select(-X) %>%
    column_to_rownames("OTUs") %>%
    select(any_of(union(meta_samples,tax_cols)))

sample_cols <- setdiff(colnames(mp4_OTU_gtdb_df), tax_cols)

mp4_OTU_gtdb_filtered <- mp4_OTU_gtdb_df %>%
  rowwise() %>%
  mutate(total_count = sum(c_across(all_of(sample_cols)))) %>%
  ungroup() %>%
  filter(total_count > 0) %>%     # removes 0
  select(-total_count)

mp4_OTU_gtdb_OTU <- mp4_OTU_gtdb_filtered %>%
    select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
    as.data.frame()

mp4_OTU_gtdb_TAX <- mp4_OTU_gtdb_filtered %>%
    select(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    as.data.frame()

mp4_OTU_gtdb_meco <- microtable$new( otu_table = mp4_OTU_gtdb_OTU, tax_table = mp4_OTU_gtdb_TAX, sample_table = metadata_df)

### SingleM OTU gtdb
singlem_OTU_gtdb_df <- read.csv(singlem_OTU_gtdb, header = TRUE, stringsAsFactors = FALSE) %>%
    select(-X) %>%
    column_to_rownames("OTUs") %>%
    select(any_of(union(meta_samples,tax_cols)))

sample_cols <- setdiff(colnames(singlem_OTU_gtdb_df), tax_cols)

singlem_OTU_gtdb_filtered <- singlem_OTU_gtdb_df %>%
  rowwise() %>%
  mutate(total_count = sum(c_across(all_of(sample_cols)))) %>%
  ungroup() %>%
  filter(total_count > 2) %>%     # removes 0-2
  filter(!grepl("Eukaryota", Kingdom)) %>%
  select(-total_count)

singlem_OTU_gtdb_OTU <- singlem_OTU_gtdb_filtered %>%
    select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
    as.data.frame()

singlem_OTU_gtdb_TAX <- singlem_OTU_gtdb_filtered %>%
    select(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    as.data.frame()

singlem_OTU_gtdb_meco <- microtable$new( otu_table = singlem_OTU_gtdb_OTU, tax_table = singlem_OTU_gtdb_TAX, sample_table = metadata_df)

# Comparison of different methods
vsearch_OTU_silva_meco
dada2_ASV_silva_meco
vsearch_gtdb_meco
dada2_ASV_gtdb_meco
mp4_OTU_sgb_meco
mp4_OTU_gtdb_meco
singlem_OTU_gtdb_meco

method_list <- c("vsearch_OTU_silva_meco","dada2_ASV_silva_meco","vsearch_OTU_gtdb_meco","dada2_ASV_gtdb_meco","mp4_OTU_sgb_meco","mp4_OTU_gtdb_meco","singlem_OTU_gtdb_meco")

## Get the summary of OTU predicted per method
### Summary OTU/ASV per method
method <- method_list[1]

method_df <- tibble(raw = method_list) %>%
  separate(
    raw, into = c("Method", "Output", "Database", "drop"),
    sep = "_", remove = FALSE
  ) %>%
  select(-drop) %>%
  mutate(
    Method = str_to_sentence(Method),
    Database = str_to_upper(Database),
    Method = recode(Method, "Mp4" = "Metaphlan4")
  )

meco_list <- setNames(mget(method_list, envir = .GlobalEnv), method_list)

#### Summary of predicted per method
total_otu_summary <- imap_dfr(
  meco_list,
  ~ tibble(
      ID = .y,
      Total_OTUs = nrow(.x$otu_table)
    )
)
total_otu_summary

#### Per sample OTU Count
otu_per_sample_tbl <- imap_dfr(
  meco_list,
  ~ {
    mat <- as.matrix(.x$otu_table)
    
    # ensure OTUs x samples
    if (nrow(mat) < ncol(mat)) mat <- t(mat)
    
    tibble(
      sample = colnames(mat),
      n_OTUs = colSums(mat > 0),
      ID = .y
    )
  }
)

otu_per_sample_tbl

#### Summarise
output_otu_summary <- otu_per_sample_tbl %>%
    pivot_wider(names_from = sample, values_from = n_OTUs) %>%
    left_join(total_otu_summary, by = c("ID" = "ID") ) %>%
    left_join(method_df, by = c("ID" = "raw")) %>%
    select(-ID) %>%                                  # drop ID
    relocate(Method, Output, Database, Total_OTUs) %>%  # force column order
    select(
        Method, Output, Database, Total_OTUs,
        sort(setdiff(names(.), c("Method","Output","Database","Total_OTUs")))
    ) %>%
    arrange(Database,Method)


#### Visualization

otu_long <- output_otu_summary %>%
  pivot_longer(
    cols = -c(Method, Output, Database, Total_OTUs),
    names_to = "Sample",
    values_to = "Observed_OTUs"
  )

plot_summary <- ggplot(otu_long, aes(Method, Observed_OTUs, fill = Database)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.5) +
  scale_y_log10() +
  theme_bw() +
  labs(
    y = "Observed OTUs / ASVs per sample (log10)",
    x = NULL,
    title = "Per-sample richness across methods"
  )

ggsave(plot = plot_summary, filename = "otu_richness_per_sample.png", width = 8, height = 6, dpi = 300)

### Get species per sample per method
#### Visualize using upset plot how many taxas are shared per method
method_list_gtdb <- method_df %>%
  filter(Database == "GTDB") %>%
  pull(raw)

meco_list_gtdb <- setNames(mget(method_list_gtdb, envir = .GlobalEnv), method_list_gtdb)

#### Extract set of lineages per method
make_lineage <- function(tx, ranks = c("Kingdom","Phylum","Class","Order","Family","Genus", "Species"),
                         sep = "; ") {
  ranks <- ranks[ranks %in% names(tx)]
  if (length(ranks) == 0) stop("No requested ranks found in tax_table.")

  tx %>%
    transmute(lineage = pmap_chr(across(all_of(ranks)), \(...) {
      parts <- c(...)
      parts <- parts[!is.na(parts) & parts != ""]
      paste(parts, collapse = sep)
    })) %>%
    mutate(lineage = na_if(lineage, "")) %>%
    filter(!is.na(lineage)) %>%
    pull(lineage)
}

taxa_sets <- imap(meco_list_gtdb, ~ {
  tx <- as_tibble(.x$tax_table, rownames = "feature_id")
  make_lineage(tx) %>%
    unique()
})

names(taxa_sets) <- c(
  "VSEARCH(GTDB)",
  "DADA2(GTDB)",
  "MetaPhlAn4(GTDB)",
  "SingleM(GTDB)"
)

#### Visualize using UpSet Plot
library(ComplexUpset)
#library(ggbreak)

all_taxa <- sort(unique(unlist(taxa_sets)))
upset_df <- tibble(taxon = all_taxa)

for (nm in names(taxa_sets)) {
  upset_df[[nm]] <- upset_df$taxon %in% taxa_sets[[nm]]
}

p_upset <- upset(
  upset_df,
  intersect = names(taxa_sets),
  name = "GTDB methods",
  base_annotations = list(
    "Species" = intersection_size(
      counts = TRUE,
      bar_number_threshold = 1,
      text = list(
        vjust = -0.4,
        size = 3
        ) 
    ) +
    theme(
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 15)
        )
  ),
  set_sizes = upset_set_size() + 
    theme(axis.text.x = element_text(size = 12, angle = 45))
)+
theme(
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5)
)

p_upset
ggsave(plot = p_upset, filename = "upset_plot_species.png", width = 8, height = 6, dpi = 300)
