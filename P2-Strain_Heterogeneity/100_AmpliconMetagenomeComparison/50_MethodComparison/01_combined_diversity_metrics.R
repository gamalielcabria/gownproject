# Create diversity
library(tidyverse)
library(microeco)

# Paths of combined otu tax files

vsearch_OTU_silva <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/20_qiime2_otu/silva_bact_arch_combined_qiimeOTU.csv"
vsearch_OTU_gtdb <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/20_qiime2_otu/gtdb_bact_arch_combined_qiimeOTU.csv"
dada2_ASV_gtdb <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv/gtdb_bact_arch_combined_qiimeASV.csv"
dada2_ASV_silva <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv/silva_bact_arch_combined_qiimeASV.csv"
mp4_OTU_sgb <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/40_metaphlan4_otu/mp4-species_bact_arch_combined_otu.csv"
singlem_OTU_gtdb <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/11_SingleM_genes/gtdb_bact_arch_combined_singlem-sum-OTU.csv"

metadata_path <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/10_GenusAggregate_Sum/sample_table_Meta.csv"

# Read files and separate into otu and tax tables
## Read and Clean metadata
metadata_df <- read.csv(metadata_path, stringsAsFactors = FALSE) %>%
    separate(X, sep = "-", into = c("Sample_Code", "Treatment", "GOWN_Well"), extra = "drop") %>%
    select(Metagenome_Code, GOWN_Well, Year, Sample_Code) %>%
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

vsearch_silva_meco <- microtable$new( otu_table = vsearch_silva_OTU, tax_table = vsearch_silva_TAX, sample_table = metadata_df)


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

vsearch_gtdb_meco <- microtable$new( otu_table = vsearch_gtdb_OTU, tax_table = vsearch_gtdb_TAX, sample_table = metadata_df)

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

### SingleM OTU gtdb