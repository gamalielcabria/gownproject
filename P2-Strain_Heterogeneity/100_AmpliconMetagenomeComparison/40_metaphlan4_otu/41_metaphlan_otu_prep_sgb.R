# Script for preparing the qiime taxonomy and otu table inputs
library(tidyverse)
library(microeco)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

# Set Working Directory~
setwd("/home/glbcabria/Workbench")

# Import OTU Counts, taxa, and metadata
## File Path of inputs
mp4_otu_tax_gtdb <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/40_metaphlan4_otu/abundance_table.mp4_merged_gtdb_relab.tsv.gz"
mp4_otu_tax_sgb <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/40_metaphlan4_otu/abundance_table.mp4_merged_abs.tsv.gz"

outpath <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/40_metaphlan4_otu"

## Importing inputs
otu_tax_sgb_mdf <- read_tsv(mp4_otu_tax_sgb, skip = 1) %>%
    rename(otu_id = `clade_name`)

unassigned_otu_tax_sgb <- otu_tax_sgb_mdf %>%
    filter(otu_id == "UNCLASSIFIED")

assigned_otu_tax_sgb <- otu_tax_sgb_mdf %>%
    filter(otu_id != "UNCLASSIFIED")


## Filter at different taxa levels
tax_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Strain")

### Selecting Species Level
# sgb_strain_otu_tax <- assigned_otu_tax_sgb %>%
#     filter(grepl("t__", otu_id)) %>%
#     separate(otu_id, into = tax_cols, sep = ";", fill = "right", extra = "drop") %>%
#     rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
#         .cols = -c(Kingdom,Phylum,Class,Order,Family,Genus,Species,Strain)
#     ) %>%
#     select(
#         all_of(tax_cols),
#         sort(setdiff(colnames(.),tax_cols))
#     ) %>%
#     mutate(OTUs = ifelse(grepl("k__Bacteria", Kingdom), paste0("bOTU_", row_number()), paste0("aOTU_", row_number())))

sgb_species_otu_tax <- assigned_otu_tax_sgb %>%
    filter(grepl("s__", otu_id)) %>%
    separate(otu_id, into = tax_cols[1:7], sep = ";", fill = "right", extra = "drop") %>%
    rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
        .cols = -c(Kingdom,Phylum,Class,Order,Family,Genus,Species)
    ) %>%
    select(
        all_of(tax_cols[1:7]),
        sort(setdiff(colnames(.),tax_cols))
    ) %>%
    mutate(OTUs = ifelse(grepl("k__Bacteria", Kingdom), paste0("bOTU_", row_number()), paste0("aOTU_", row_number()))) %>%
    bind_rows(
        unassigned_otu_tax_sgb %>%
        separate(otu_id, into = tax_cols[1:7], sep = ";", fill = "right", extra = "drop") %>%
        rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
            .cols = -c(Kingdom,Phylum,Class,Order,Family,Genus,Species)
        ) %>%
        select(
            all_of(tax_cols[1:7]),
            sort(setdiff(colnames(.),tax_cols))
        ) %>%
        mutate(OTUs = paste0("uOTU_", row_number()))
    )

sgb_genus_otu_tax <- assigned_otu_tax_sgb %>%
    filter(grepl("g__", otu_id) & !grepl("s__", otu_id)) %>%
    separate(otu_id, into = tax_cols[1:6], sep = ";", fill = "right", extra = "drop") %>%
    rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
        .cols = -c(Kingdom,Phylum,Class,Order,Family,Genus)
    ) %>%
    select(
        all_of(tax_cols[1:6]),
        sort(setdiff(colnames(.),tax_cols))
    ) %>%
    mutate(OTUs = ifelse(grepl("k__Bacteria", Kingdom), paste0("bOTU_", row_number()), paste0("aOTU_", row_number()))) %>%
    bind_rows(
        unassigned_otu_tax_sgb %>%
        separate(otu_id, into = tax_cols[1:6], sep = ";", fill = "right", extra = "drop") %>%
        rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
            .cols = -c(Kingdom,Phylum,Class,Order,Family,Genus,Species)
        ) %>%
        select(
            all_of(tax_cols[1:6]),
            sort(setdiff(colnames(.),tax_cols))
        ) %>%
        mutate(OTUs = paste0("uOTU_", row_number()))
    )

sgb_family_otu_tax <- assigned_otu_tax_sgb %>%
    filter(grepl("f__", otu_id) & !grepl("g__", otu_id) & !grepl("s__", otu_id)) %>%
    separate(otu_id, into = tax_cols[1:5], sep = ";", fill = "right", extra = "drop") %>%
    rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
        .cols = -c(Kingdom,Phylum,Class,Order,Family)
    ) %>%
    select(
        all_of(tax_cols[1:5]),
        sort(setdiff(colnames(.),tax_cols))
    ) %>%
    mutate(OTUs = ifelse(grepl("k__Bacteria", Kingdom), paste0("bOTU_", row_number()), paste0("aOTU_", row_number()))) %>%
    bind_rows(
        unassigned_otu_tax_sgb %>%
        separate(otu_id, into = tax_cols[1:5], sep = ";", fill = "right", extra = "drop") %>%
        rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
            .cols = -c(Kingdom,Phylum,Class,Order,Family,Genus,Species)
        ) %>%
        select(
            all_of(tax_cols[1:5]),
            sort(setdiff(colnames(.),tax_cols))
        ) %>%
        mutate(OTUs = paste0("uOTU_", row_number()))
    )

sgb_phylum_otu_tax <- assigned_otu_tax_sgb %>%
    filter(grepl("p__", otu_id) & !grepl("c__", otu_id)) %>%
    separate(otu_id, into = tax_cols[1:2], sep = ";", fill = "right", extra = "drop") %>%
    rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
        .cols = -c(Kingdom,Phylum)
    ) %>%
    select(
        all_of(tax_cols[1:2]),
        sort(setdiff(colnames(.),tax_cols))
    ) %>%
    mutate(OTUs = ifelse(grepl("k__Bacteria", Kingdom), paste0("bOTU_", row_number()), paste0("aOTU_", row_number()))) %>%
    bind_rows(
        unassigned_otu_tax_sgb %>%
        separate(otu_id, into = tax_cols[1:2], sep = ";", fill = "right", extra = "drop") %>%
        rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
            .cols = -c(Kingdom,Phylum,Class,Order,Family,Genus,Species)
        ) %>%
        select(
            all_of(tax_cols[1:2]),
            sort(setdiff(colnames(.),tax_cols))
        ) %>%
        mutate(OTUs = paste0("uOTU_", row_number()))
    )

### Fill missing Taxa
### Fill in missing taxonomic ranks
tax_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
filled_sgb_species_otu_tax <- sgb_species_otu_tax %>%
    mutate(Kingdom = ifelse(Kingdom == "Unassigned", "d__Unassigned", Kingdom)) %>%
    mutate(
        across(
        Kingdom:Species,
        ~ {
            x <- replace_na(.x, "")
            core <- str_remove(x, "^[a-z]__")
            if_else(core %in% c("", "unclassified"), "unclassified", x)
        }
        )
    )

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

filled_sgb_species_otu_tax2 <- filled_sgb_species_otu_tax
filled_sgb_species_otu_tax2[tax_cols] <- t(apply(filled_sgb_species_otu_tax2[tax_cols],1, fill_from_parent))


## Save as csv files
write.csv(filled_sgb_species_otu_tax2, file = paste0(outpath, "/mp4-species_bact_arch_combined_otu.csv"))
# write.csv(sgb_genus_otu_tax, file = paste0(outpath, "/mp4-genus_bact_arch_combined_otu.csv"))
# write.csv(sgb_family_otu_tax, file = paste0(outpath, "/mp4-family_bact_arch_combined_otu.csv"))
# write.csv(sgb_phylum_otu_tax, file = paste0(outpath, "/mp4-phylum_bact_arch_combined_otu.csv"))



###################
# PROCESISNG GTDB #
###################
## Importing inputs
otu_tax_gtdb_mdf <- read_tsv(mp4_otu_tax_gtdb, skip = 1) %>%
    rename(otu_id = `clade_name`)

unassigned_otu_tax_gtdb <- otu_tax_gtdb_mdf %>%
    filter(otu_id == "UNCLASSIFIED")

assigned_otu_tax_gtdb <- otu_tax_gtdb_mdf %>%
    filter(otu_id != "UNCLASSIFIED")


## Filter at different taxa levels
tax_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Strain")

gtdb_species_otu_tax <- assigned_otu_tax_gtdb %>%
    filter(grepl("s__", otu_id)) %>%
    separate(otu_id, into = tax_cols[1:7], sep = ";", fill = "right", extra = "drop") %>%
    rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
        .cols = -c(Kingdom,Phylum,Class,Order,Family,Genus,Species)
    ) %>%
    select(
        all_of(tax_cols[1:7]),
        sort(setdiff(colnames(.),tax_cols))
    ) %>%
    mutate(OTUs = ifelse(grepl("k__Bacteria", Kingdom), paste0("bOTU_", row_number()), paste0("aOTU_", row_number()))) %>%
    bind_rows(
        unassigned_otu_tax_sgb %>%
        separate(otu_id, into = tax_cols[1:7], sep = ";", fill = "right", extra = "drop") %>%
        rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
            .cols = -c(Kingdom,Phylum,Class,Order,Family,Genus,Species)
        ) %>%
        select(
            all_of(tax_cols[1:7]),
            sort(setdiff(colnames(.),tax_cols))
        ) %>%
        mutate(OTUs = paste0("uOTU_", row_number()))
    )

gtdb_genus_otu_tax <- assigned_otu_tax_gtdb %>%
    filter(grepl("g__", otu_id) & !grepl("s__", otu_id)) %>%
    separate(otu_id, into = tax_cols[1:6], sep = ";", fill = "right", extra = "drop") %>%
    rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
        .cols = -c(Kingdom,Phylum,Class,Order,Family,Genus)
    ) %>%
    select(
        all_of(tax_cols[1:6]),
        sort(setdiff(colnames(.),tax_cols))
    ) %>%
    mutate(OTUs = ifelse(grepl("k__Bacteria", Kingdom), paste0("bOTU_", row_number()), paste0("aOTU_", row_number()))) %>%
    bind_rows(
        unassigned_otu_tax_gtdb %>%
        separate(otu_id, into = tax_cols[1:6], sep = ";", fill = "right", extra = "drop") %>%
        rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
            .cols = -c(Kingdom,Phylum,Class,Order,Family,Genus)
        ) %>%
        select(
            all_of(tax_cols[1:6]),
            sort(setdiff(colnames(.),tax_cols))
        ) %>%
        mutate(OTUs = paste0("uOTU_", row_number()))
    )

gtdb_family_otu_tax <- assigned_otu_tax_gtdb %>%
    filter(grepl("f__", otu_id) & !grepl("g__", otu_id) & !grepl("s__", otu_id)) %>%
    separate(otu_id, into = tax_cols[1:5], sep = ";", fill = "right", extra = "drop") %>%
    rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
        .cols = -c(Kingdom,Phylum,Class,Order,Family)
    ) %>%
    select(
        all_of(tax_cols[1:5]),
        sort(setdiff(colnames(.),tax_cols))
    ) %>%
    mutate(OTUs = ifelse(grepl("k__Bacteria", Kingdom), paste0("bOTU_", row_number()), paste0("aOTU_", row_number()))) %>%
    bind_rows(
        unassigned_otu_tax_gtdb %>%
        separate(otu_id, into = tax_cols[1:5], sep = ";", fill = "right", extra = "drop") %>%
        rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
            .cols = -c(Kingdom,Phylum,Class,Order,Family)
        ) %>%
        select(
            all_of(tax_cols[1:5]),
            sort(setdiff(colnames(.),tax_cols))
        ) %>%
        mutate(OTUs = paste0("uOTU_", row_number()))
    )

gtdb_phylum_otu_tax <- assigned_otu_tax_gtdb %>%
    filter(grepl("p__", otu_id) & !grepl("c__", otu_id)) %>%
    separate(otu_id, into = tax_cols[1:2], sep = ";", fill = "right", extra = "drop") %>%
    rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
        .cols = -c(Kingdom,Phylum)
    ) %>%
    select(
        all_of(tax_cols[1:2]),
        sort(setdiff(colnames(.),tax_cols))
    ) %>%
    mutate(OTUs = ifelse(grepl("k__Bacteria", Kingdom), paste0("bOTU_", row_number()), paste0("aOTU_", row_number()))) %>%
    bind_rows(
        unassigned_otu_tax_gtdb %>%
        separate(otu_id, into = tax_cols[1:2], sep = ";", fill = "right", extra = "drop") %>%
        rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
            .cols = -c(Kingdom,Phylum)
        ) %>%
        select(
            all_of(tax_cols[1:2]),
            sort(setdiff(colnames(.),tax_cols))
        ) %>%
        mutate(OTUs = paste0("uOTU_", row_number()))
    )

### Fill missing Taxa
### Fill in missing taxonomic ranks
tax_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
filled_gtdb_species_otu_tax <- gtdb_species_otu_tax %>%
    mutate(Kingdom = ifelse(Kingdom == "Unassigned", "d__Unassigned", Kingdom)) %>%
    mutate(
        across(
        Kingdom:Species,
        ~ {
            x <- replace_na(.x, "")
            core <- str_remove(x, "^[a-z]__")
            if_else(core %in% c("", "unclassified"), "unclassified", x)
        }
        )
    )

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

filled_gtdb_species_otu_tax2 <- filled_gtdb_species_otu_tax
filled_gtdb_species_otu_tax2[tax_cols] <- t(apply(filled_gtdb_species_otu_tax2[tax_cols],1, fill_from_parent))

write.csv(filled_gtdb_species_otu_tax2, file = paste0(outpath, "/mp4-species_bact_arch_combined_relabgtdb.csv"))
# write.csv(gtdb_genus_otu_tax, file = paste0(outpath, "/mp4-genus_bact_arch_combined_relabgtdb.csv"))
# write.csv(gtdb_family_otu_tax, file = paste0(outpath, "/mp4-family_bact_arch_combined_relabgtdb.csv"))
# write.csv(gtdb_phylum_otu_tax, file = paste0(outpath, "/mp4-phylum_bact_arch_combined_relabgtdb.csv"))

#################
# END OF SCRIPT #
#################
