# Create diversity
library(tidyverse)
library(microeco)

# Paths of combined otu tax files

vsearch_OTU_silva <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/20_qiime2_otu/silva_bact_arch_combined_qiimeOTU.csv"
vsearch_OTU_gtdb <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/20_qiime2_otu/gtdb_bact_arch_combined_qiimeOTU.csv"
dada2_ASV_gtdb <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv/gtdb_bact_arch_combined_qiimeASV.csv"
dada2_ASV_silva <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/30_qiime2_asv/silva_bact_arch_combined_qiimeASV.csv"
mp4_OTU_sgb <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/40_metaphlan4_otu/mp4-species_bact_arch_combined_otu.csv"
singlem_OTU_genus <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/11_SingleM_genes/gtdb_bact_arch_combined_singlem-sum-OTU.csv"

metadata_path <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/10_GenusAggregate_Sum/sample_table_Meta.csv"

# Read files and separate into otu and tax tables
## Read and Clean metadata
metadata_df <- read.csv(metadata_path, stringsAsFactors = FALSE) %>%
    separate(X, sep = "-", into = c("Sample_Code", "Treatment", "GOWN_Well")) %>%
    select(Metagenome_Code, GOWN_Well, Year, Sample_Code)

## Read and Break into component taxa
vsearch_OTU_silva_df <- read.csv(vsearch_OTU_silva, stringsAsFactors = FALSE) %>%
    select(-X) %>%
    column_to_rownames("OTUs")

vsearch_silva_OTU <- vsearch_OTU_silva_df %>%
    select()