library("tidyverse")
library("microeco")

# Import OTU Counts, taxa, and metadata
## File Path of inputs
amp_combined_taxa_otu <- "/home/glbcabria/Workbench/gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/10_GenusAggregate_Sum/otu_taxa_combined_Amp.csv"
mtg_combined_taxa_otu <- "/home/glbcabria/Workbench/gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/10_GenusAggregate_Sum/otu_taxa_combined_Meta.csv"
amp_sample_data <- "/home/glbcabria/Workbench/gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/10_GenusAggregate_Sum/sample_table_Amp_edited.csv"
mtg_sample_data <- "/home/glbcabria/Workbench/gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/10_GenusAggregate_Sum/sample_table_Meta.csv"

## Importing inputs
combined_amp <- read.csv(amp_combined_taxa_otu) %>% rename(ASVs=X)
combined_mtg <- read.csv(mtg_combined_taxa_otu) %>% rename(ASVs=X)

mtg_metadata <- read.csv(mtg_sample_data) %>% rename(MtgCols = X) %>%
    separate(MtgCols, into = c("Extra1","Extra2","GOWN_Well", "Extra4"), remove = FALSE) %>%
    select(-Extra1, -Extra2, -Extra4)
amp_metadata <- read.csv(amp_sample_data) %>% rename(AmpCols = X) %>%
    left_join(mtg_metadata %>% select(Metagenome_Code, GOWN_Well), by = "Metagenome_Code") 

# Process Sample Data
## Create an Amp meco object

### Split the combined_amp to bact and arch tables
combined_zero_ASVs <- combined_amp %>%
  rowwise() %>%
  filter(!(sum(c_across(DMO.GOWN.C1.20200316:MDI.GOWN.E9.20200316), na.rm = TRUE) == 0)) %>%
  ungroup()

MDI_combined_amp <- combined_amp %>%
    select(ASVs,Kingdom,Phylum,Class,Order,Family,Genus,starts_with("MDI")) %>%
    rename_with(
        ~ str_split(.x, "\\.", simplify = TRUE)[,3],
        .cols = -c(ASVs,Kingdom,Phylum,Class,Order,Family,Genus)
    ) %>%
    rowwise() %>%
    filter(!(sum(c_across(C1:E9), na.rm = TRUE) == 0) ) %>%
    ungroup()

DMO_combined_amp <- combined_amp %>%
    select(ASVs,Kingdom,Phylum,Class,Order,Family,Genus,starts_with("DMO")) %>%
    rename_with(
        ~ str_split(.x, "\\.", simplify = TRUE)[,3],
        .cols = -c(ASVs,Kingdom,Phylum,Class,Order,Family,Genus)
    )%>%
    rowwise() %>%
    filter(!(sum(c_across(C1:E9), na.rm = TRUE) == 0)) %>%
    ungroup()

# identical(
#   DMO_combined_amp %>% arrange(ASVs) %>% select(ASVs, Kingdom:Genus),
#   MDI_combined_amp %>% arrange(ASVs) %>% select(ASVs, Kingdom:Genus)
# )

rearranged_combined_amp <- bind_rows(DMO_combined_amp, MDI_combined_amp) %>%
    group_by(ASVs, Kingdom, Phylum, Class, Order, Family, Genus) %>%
    summarise(
        across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
        .groups = "drop"
    )
###
taxa_table_amp <- combined_amp %>%
    select(ASVs,Kingdom,Phylum,Class,Order,Family,Genus) %>%
    column_to_rownames("ASVs")

otu_table_amp<- combined_amp %>%
    select(-Kingdom,-Phylum,-Class,-Order,-Family,-Genus) 
