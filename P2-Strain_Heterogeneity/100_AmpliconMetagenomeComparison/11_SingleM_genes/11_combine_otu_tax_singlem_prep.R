library(tidyverse)

# inputs
df <- read_tsv("gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/11_SingleM_genes/combined_singlem_otu_table.tsv")

outpath <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/11_SingleM_genes/"

# Sum or Mean the counts
summed_df <- df %>%
  group_by(sample, taxonomy) %>%
  summarise(
    num_hits = sum(num_hits, na.rm = TRUE),
    .groups = "drop"
  )

otu_wide <- summed_df %>%
    pivot_wider(
    names_from  = sample,
    values_from = num_hits,
    values_fill = 0
    ) %>%
    mutate(
    OTUs = case_when(
      str_detect(taxonomy, regex("bacteria", ignore_case = TRUE)) ~
        paste0("bOTU_", row_number()),
      str_detect(taxonomy, regex("archaea", ignore_case = TRUE)) ~
        paste0("aOTU_", row_number()),
      TRUE ~ paste0("uOTU_", row_number())
        )
    ) %>%
    separate(
        taxonomy,
        into = c("root","Kingdom","Phylum","Class","Order","Family","Genus","Species"),
        sep = "; ",
        fill = "right",
        remove = FALSE
    ) %>%
    mutate(
    Kingdom = if_else(taxonomy == "Root", "d__Unassigned", Kingdom)
    ) %>%
    select(-taxonomy, -root) %>%
    rename_with( ~ str_split(.x, "-", simplify = TRUE)[,2], .cols = -c(OTUs,Kingdom,Phylum,Class,Order,Family,Genus,Species)    )

## Separate OTU and TAX for 
singlem_taxa <- otu_wide %>%
    select(OTUs, Kingdom, Phylum, Class, Order, Family, Genus, Species)

singlem_otu <- otu_wide %>%
    select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus, -Species)


### Fill in missing taxonomic ranks
tax_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
singlem_taxa_filled <- singlem_taxa %>%
    mutate(across(-OTUs, ~ replace_na(.x,"unclassified")))

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

singlem_taxa_filled2 <- singlem_taxa_filled
singlem_taxa_filled2[tax_cols] <- t(apply(singlem_taxa_filled2[tax_cols],1, fill_from_parent))

# Recombine OTU and Tax tables

combined_otu_tax_singlem <- singlem_taxa_filled2 %>%
    left_join(singlem_otu, by = "OTUs")

# Save the file to csv
# write.csv(combined_otu_tax_singlem, file = paste0(outpath, "gtdb_bact_arch_combined_singlem-sum-OTU.csv"))