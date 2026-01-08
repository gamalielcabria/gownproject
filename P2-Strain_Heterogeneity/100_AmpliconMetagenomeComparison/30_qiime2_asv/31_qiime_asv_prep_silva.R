# Script for preparing the qiime taxonomy and otu table inputs
library(tidyverse)
library(microeco)

# Set Working Directory~
setwd("/home/glbcabria/Workbench")

# Import OTU Counts, taxa, and metadata
## File Path of inputs
qiime_otu_silva <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv/01_qiime2_asv_silva/asv-table.tsv.gz"
qiime_tax_silva <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv/01_qiime2_asv_silva/taxonomy.tsv.gz"

outpath <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv"

## Importing inputs
otu_silva_qdf <- read_tsv(qiime_otu_silva, skip = 1) %>%
    rename(otu_id = `#OTU ID`)

tax_silva_qdf <- read_tsv(qiime_tax_silva, col_names = c("otu_id","taxa","fraction")) %>%
    select(-fraction) %>%
    separate(taxa, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";", fill = "right", extra = "drop") 

### Fill in missing taxonomic ranks
tax_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_silva_qdf_filled <- tax_silva_qdf %>%
    mutate(Kingdom = ifelse(Kingdom == "Unassigned", "d__Unassigned", Kingdom)) %>%
    mutate(across(-otu_id, ~ replace_na(.x,"unclassified")))

fill_from_parent <- function(row) {
  #core <- strip_prefix(row)
  is_uncl <- tolower(row) == "unclassified"

  for (i in which(is_uncl)) {
    if (i == 1) next
    j <- max(which(
      !is.na(row[1:(i-1)]) &
      row[1:(i-1)] != "" &
      tolower(row[1:(i-1)]) != "unclassified"
    ), na.rm = TRUE)

    if (is.finite(j)) row[i] <- row[j]
  }
  row
}

tax_silva_qdf_filled2 <- tax_silva_qdf_filled
tax_silva_qdf_filled2[tax_cols] <- t(apply(tax_silva_qdf_filled2[tax_cols],1, fill_from_parent))

## Combining OTU and Taxa tables for unified OTU_ID
combined_otu_tax_silva <- otu_silva_qdf %>%
    left_join(tax_silva_qdf_filled2, by = "otu_id") %>%
    mutate(OTUs = paste0("OTU_", row_number())) %>%
    rename_with(~gsub("-",".",.x)) %>%
    select(-otu_id)




###################################
# SEPARATE BACT AND ARCH ANALYSIS #
###################################

# Analysis of silva OTUs

## Get the qiime OTU (qOTU) for Bacteria using the MDI tag
tax_cols <- c("OTUs", "Kingdom", "Phylum", "Class", "Order",
              "Family", "Genus", "Species")

bact_combined_qOTU <- combined_otu_tax_silva %>%
    select(OTUs,Kingdom,Phylum,Class,Order,Family,Genus,Species,starts_with("MDI")) %>%
    rename_with(~ str_split(.x, "\\.", simplify = TRUE)[,4],
        .cols = -c(OTUs,Kingdom,Phylum,Class,Order,Family,Genus,Species)
    ) %>%
    select(
        all_of(tax_cols),
        sort(setdiff(colnames(.),tax_cols))
    ) %>%
    filter( grepl("Bacteria|Unassigned", Kingdom) ) %>%
    rowwise() %>%
    filter(!(sum(c_across(C1:E9), na.rm = TRUE) == 0) ) %>%
    ungroup() %>%
    mutate(OTUs = gsub("OTU_", "bOTU_", OTUs))

arch_combined_qOTU <- combined_otu_tax_silva %>%
    select(OTUs,Kingdom,Phylum,Class,Order,Family,Genus,Species,starts_with("DMO")) %>%
    rename_with(~ str_split(.x, "\\.", simplify = TRUE)[,4],
        .cols = -c(OTUs,Kingdom,Phylum,Class,Order,Family,Genus,Species)
    ) %>%
    select(
        all_of(tax_cols),
        sort(setdiff(colnames(.),tax_cols))
    ) %>%
    filter( grepl("Archaea|Unassigned", Kingdom) ) %>%
    rowwise() %>%
    filter(!(sum(c_across(C1:E9), na.rm = TRUE) == 0) ) %>%
    ungroup() %>%
    mutate(OTUs = gsub("OTU_", "aOTU_", OTUs))

# unassigned_combined_qOTU <- combined_otu_tax_silva %>%
#     filter( grepl("Unassigned", Kingdom) ) 
#     nrow(unassigned_combined_qOTU)


silva_bact_arch_combined_ASV <- bind_rows(bact_combined_qOTU, arch_combined_qOTU)

## Write the combined OTU and Taxa tables to CSV files
write.csv(silva_bact_arch_combined_ASV, file = paste0(outpath, "/silva_bact_arch_combined_qiimeASV.csv")) 
# silva_bact_arch_combined_OTU2 <- read.csv(paste0(outpath, "/silva_bact_arch_combined_qiimeOTU.csv")) %>% select(-X)
# bact_combined_qOTU2 <- silva_bact_arch_combined_OTU2 %>% filter( grepl("bOTU", OTUs) )
# arch_combined_qOTU2 <- silva_bact_arch_combined_OTU2 %>% filter( grepl("aOTU", OTUs) )

# ## Separate OTU and Taxa tables
# otu_table_silva_all <- silva_bact_arch_combined_OTU2 %>%
#     select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
#     column_to_rownames(var = "OTUs")
# tax_table_silva_all <- silva_bact_arch_combined_OTU %>%
#     select(OTUs, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
#     column_to_rownames(var = "OTUs")

# otu_table_silva_bac <- bact_combined_qOTU2 %>%
#     select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
#     column_to_rownames(var = "OTUs")
# tax_table_silva_bac <- bact_combined_qOTU2 %>%
#     select(OTUs, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
#     column_to_rownames(var = "OTUs")

# otu_table_silva_arch <- arch_combined_qOTU2 %>%
#     select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species) %>%
#     column_to_rownames(var = "OTUs")
# tax_table_silva_arch <- arch_combined_qOTU2 %>%
#     select(OTUs, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
#     column_to_rownames(var = "OTUs")
