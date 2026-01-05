# Script for preparing the qiime taxonomy and otu table inputs
library(tidyverse)
library(microeco)

# Set Working Directory
#setwd("/home/glbcabria/Workbench")

# Import OTU Counts, taxa, and metadata
## File Path of inputs
qiime_otu_gtdb <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/20_qiime2_otu/02_qiime2_otu_gtdb/otu-table.tsv.gz"
qiime_otu_silva <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/20_qiime2_otu/01_qiime2_otu_silva/otu-table.tsv.gz"
qiime_tax_gtdb <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/20_qiime2_otu/02_qiime2_otu_gtdb/taxonomy.tsv.gz"
qiime_tax_silva <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/20_qiime2_otu/01_qiime2_otu_silva/taxonomy.tsv.gz"

outpath <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/20_qiime2_otu"

## Importing inputs
otu_gtdb_qdf <- read_tsv(qiime_otu_gtdb, skip = 1) %>%
    rename(otu_id = `#OTU ID`)

tax_gtdb_qdf <- read_tsv(qiime_tax_gtdb, col_names = c("otu_id","taxa","fraction")) #%>%
    select(-fraction)

combined_otu_tax_gtdb <- otu_gtdb_qdf %>%
    left_join(tax_gtdb_qdf, by = "otu_id") %>%
    separate(taxa, into = c("kingdom","phylum","class","order","family","genus","species"), sep = "; ", fill = "right", extra = "drop") %>%
    mutate(OTU = paste0("OTU_", row_number())) %>%
    rename_with(~gsub("-",".",.x)) %>%
    select(-otu_id)

# combined_zero_OTUs_gtdb <- combined_otu_tax_gtdb %>%
#   rowwise() %>%
#   filter(!(sum(c_across(DMO.arch.GOWN.C1.20200316:MDI.bact.GOWN.C5.20200316), na.rm = TRUE) == 0)) %>%
#   ungroup()
