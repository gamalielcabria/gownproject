# Processing combined OTU TAXA tables into Microeco Objects
library(tidyverse)
library(microeco)

# Setting up Work Directory (For non-combined github repo)
setwd("/home/glbcabria/Workbench/")

# DECLARING INPUTS
combined_asv_taxa <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv/gtdb_bact_arch_combined_qiimeASV_GOWN24-25.csv"
metadata_fp <- "gownproject/P0-GOWN/200_denitrificationpotential/00_rawdata/run251210_Metadata_amplicon.new.csv"
#metadata_fp <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv/gtdb_bact_arch_combined_qiimeASV_GOWN24-25.metadata.csv"

# READING INPUTS
dat <- readr::read_csv(combined_asv_taxa, show_col_types = FALSE)
meta <- readr::read_csv(metadata_fp, show_col_types = FALSE)

## Drop the leading index column if present (often named "" or "...1")
idx_cols <- c("", "...1")
dat <- dat %>% select(-any_of(idx_cols))

meta <- meta %>%
  select(-any_of(idx_cols)) %>%
  column_to_rownames("Filename")

## Safety Checks
### Do Taxa data included?
stopifnot("OTUs" %in% colnames(dat))

tax_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
missing_tax <- setdiff(tax_cols, colnames(dat))
if (length(missing_tax) > 0) {
  stop("Missing expected taxonomy columns: ", paste(missing_tax, collapse = ", "))
}

### Do sample columnds included?
# --- 2) Identify sample columns (everything after Species) ---
species_pos <- which(colnames(dat) == "Species")
sample_cols <- colnames(dat)[(species_pos + 1):ncol(dat)]
if (length(sample_cols) == 0) stop("No sample columns detected after 'Species'.")


# BUILDING OTU AND TAXA TABLES
## Removing Singletons and Doubletons
otu_tmp <- dat %>%
  select(OTUs, all_of(sample_cols)) %>%
  column_to_rownames("OTUs") %>%
  as.matrix()
mode(otu_tmp) <- "numeric"

max_abund <- apply(otu_tmp, 1, max, na.rm = TRUE)
keep <- max_abund > 2
keep[is.na(keep)] <- FALSE
otu_filtered <- otu_tmp[keep, , drop = FALSE]

## Changing G1 to G01 and maintaining M7E the same and similar changes
colnames(otu_filtered) <- sub("^([A-Z])(\\d)$", "\\10\\2", colnames(otu_filtered))

## Build a Tax Table
tax_df <- dat %>%
  select(OTUs, all_of(tax_cols)) %>%
  distinct(OTUs, .keep_all = TRUE) %>%
  column_to_rownames("OTUs")

### Optional: clean GTDB prefixes like d__/p__/c__/...
## strip_prefix <- function(x) sub("^[a-z]__+", "", x)
## tax_df <- tax_df %>% mutate(across(all_of(tax_cols), ~ ifelse(is.na(.x), NA_character_, strip_prefix(.x))))

tax_filtered <- tax_df[rownames(otu_filtered), , drop = FALSE]
tax_mat <- as.matrix(tax_filtered)

###############################
## Building a Metadata Table ##
###############################

### ~~ SAFETY CHECK ~~ Measure how many samples only have BACT and ARCH and BOTH
#### Now proceed as before: categorize by presence in BACT vs ARCH
presence <- meta %>% 
  rename(row_id = GenomeCode) %>%
  distinct(row_id, Domain) %>%          # row_id x Domain unique combos
  mutate(flag = TRUE) %>%
  pivot_wider(names_from = Domain, values_from = flag, values_fill = FALSE) %>%
  mutate(Category = case_when(
    BACT & ARCH ~ "BOTH",
    BACT & !ARCH ~ "BACT only",
    !BACT & ARCH ~ "ARCH only",
    TRUE ~ "Other"   # e.g., only PROK or other domains
  ))

#### Attach categories back to meta
meta_categorized <- meta %>% rename(row_id = GenomeCode) %>%
  inner_join(presence %>% select(row_id, Category), by = "row_id")

#### Counts per category (only the three of interest)
category_counts <- presence %>%
  filter(Category %in% c("BACT only", "ARCH only", "BOTH")) %>%
  count(Category, name = "n") %>%
  arrange(match(Category, c("BACT only", "ARCH only", "BOTH")))

#### Lists
# bact_only_list <- meta_categorized %>% filter(Category == "BACT only") %>%
#   distinct(row_id, Well, Year) %>% arrange(row_id, Well, Year)
# arch_only_list <- meta_categorized %>% filter(Category == "ARCH only") %>%
#   distinct(row_id, Well, Year) %>% arrange(row_id, Well, Year)
# both_list <- meta_categorized %>% filter(Category == "BOTH") %>%
#   distinct(row_id, Domain, Well, Year) %>% arrange(row_id, Domain, Well, Year)

### Create a Meta table

meta_summary <- meta %>%
  mutate(Domain = toupper(str_trim(Domain))) %>%     # normalize Domain
  group_by(GenomeCode) %>%
  summarise(
    Domain  = paste(sort(unique(Domain)), collapse = ","),
    Well    = paste(sort(unique(Well)),   collapse = ","),
    Year    = paste(sort(unique(Year)),   collapse = ","),
    Project = paste(sort(unique(Project)),collapse = ","),
    .groups = "drop"
  ) %>%
  column_to_rownames("GenomeCode")

### ~~ SANITY CHECK ~~ Matching OTU/sample columns and metadata row
otu_samples  <- colnames(otu_filtered)
meta_samples <- rownames(meta_summary)
shared <- intersect(otu_samples, meta_samples)
notshared <- setdiff(meta_samples, otu_samples)
if (length(shared) == 0) stop("No overlapping sample IDs between OTU table and metadata!")

### Matching Meta and OTU samples
otu_filtered <- otu_filtered[, shared, drop = FALSE]
meta_summary <- meta_summary[shared, , drop = FALSE]

### Setting Metadata Variable Types/Class
meta_summary <- meta_summary %>%
  as.data.frame() %>%
  mutate(
    Year = as.integer(Year),
    Experiment = as.factor(Project),
    Domain = as.factor(Domain),
    Well = as.factor(Well)#, User = as.factor(User)
  )


## BUILDING FINAL TABLES
otu_table <- as.data.frame(otu_filtered)
tax_table <- as.data.frame(tax_mat)
sample_table <- meta_summary

# CREATING MICROECO OBJECT
dada2_ASV_gtdb_meco.2024.2025 <- microtable$new(
  otu_table    = otu_table,
  tax_table    = tax_table,
  sample_table = sample_table
)

saveRDS(
  dada2_ASV_gtdb_meco.2024.2025,
  file = "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv/dada2_ASV_gtdb_meco.2024.2025.rds"
)
