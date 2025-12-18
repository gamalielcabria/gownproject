library("tidyverse")
library("microeco")
library("mecodev")

# Import OTU Counts, taxa, and metadata
## File Path of inputs
amp_combined_taxa_otu <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/10_GenusAggregate_Sum/otu_taxa_combined_Amp.csv"
mtg_combined_taxa_otu <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/10_GenusAggregate_Sum/otu_taxa_combined_Meta.csv"
amp_sample_data <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/10_GenusAggregate_Sum/sample_table_Amp_edited.csv"
mtg_sample_data <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/10_GenusAggregate_Sum/sample_table_Meta.csv"

outpath <- "gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/10_GenusAggregate_Sum"

## Importing inputs
combined_amp <- read.csv(amp_combined_taxa_otu) %>% rename(ASVs=X)
combined_mtg <- read.csv(mtg_combined_taxa_otu) %>% rename(ASVs=X)

mtg_metadata <- read.csv(mtg_sample_data) %>% rename(MtgCols = X) %>%
    separate(MtgCols, into = c("Extra1","Extra2","GOWN_Well", "Extra4"), remove = FALSE) %>%
    select(-Extra1, -Extra2, -Extra4)
amp_metadata <- read.csv(amp_sample_data) %>% rename(AmpCols = X) %>%
    left_join(mtg_metadata %>% select(Metagenome_Code, GOWN_Well), by = "Metagenome_Code") 

metadata <- amp_metadata %>%
    filter(Creator == "DMO") %>%
    select(-AmpCols, -Creator) %>%
    column_to_rownames("Metagenome_Code")

#####################################
#  Process Sample Data in two ways: #
#####################################
# Mean the counts from two primers or 
# Remove the untargeted Taxa per primer

## Create an Amp meco object for Mean of ASV Counts
### Split the combined_amp to bact and arch tables
### Remove any Zeroes 
combined_zero_ASVs <- combined_amp %>%
  rowwise() %>%
  filter(!(sum(c_across(DMO.GOWN.C1.20200316:MDI.GOWN.E9.20200316), na.rm = TRUE) == 0)) %>%
  ungroup()
### Just separate MDI and DMO and remove 0s. MDI is Bact but any Archaea with counts is still retained
MDI_combined_amp <- combined_amp %>%
    select(ASVs,Kingdom,Phylum,Class,Order,Family,Genus,starts_with("MDI")) %>%
    rename_with(
        ~ str_split(.x, "\\.", simplify = TRUE)[,3],
        .cols = -c(ASVs,Kingdom,Phylum,Class,Order,Family,Genus)
    ) %>%
    rowwise() %>%
    filter(!(sum(c_across(C1:E9), na.rm = TRUE) == 0) ) %>%
    ungroup()
### DMO is Archaea Primers but 724 out of 786 is Bacteria!
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

### Here I calculate the mean for ASVs that are common between MDI and DMO
### Bacterial ASVs that are found for both will have the counts averaged between
### the counts. But if an ASV is unique only to each dataset, it should retain its
rearranged_mean_combined_amp <- bind_rows(DMO_combined_amp, MDI_combined_amp) %>%
    group_by(ASVs, Kingdom, Phylum, Class, Order, Family, Genus) %>%
    summarise(
        across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
    )

taxa_table_amp <- rearranged_mean_combined_amp %>%
    select(ASVs,Kingdom,Phylum,Class,Order,Family,Genus) %>%
    column_to_rownames("ASVs")

otu_table_amp <- rearranged_mean_combined_amp %>%
    select(-Kingdom,-Phylum,-Class,-Order,-Family,-Genus) %>%
    column_to_rownames("ASVs")

### Create the meco object for mean counts
meco_amp_mean <- microtable$new(
    otu_table = otu_table_amp,
    tax_table = taxa_table_amp,
    sample_table = metadata
)
meco_amp_mean$cal_alphadiv(measures=c("Observed","Shannon","Simpson","Pielou"))

### Identify total counts per sample
total_counts_amp <- colSums(meco_amp_mean$otu_table)
min(total_counts_amp)
max(total_counts_amp)

meco_amp_mean_rarefy_sh <- trans_rarefy$new(meco_amp_mean, alphadiv = "Shannon",  depth = c(0, 25, 100, 400, 800, 1600, 2400, 4800, 8000, 16000))
plot_amp_mean_rarefy_shannon <- meco_amp_mean_rarefy_sh$plot_rarefy()

meco_amp_mean_rarefy_simpson <- trans_rarefy$new(meco_amp_mean, alphadiv = "Simpson",  depth = c(0, 25, 100, 400, 800, 1600, 2400, 4800, 8000, 16000))
plot_amp_mean_rarefy_simpson <- meco_amp_mean_rarefy_simpson$plot_rarefy()

meco_amp_mean_rarefy_pielou <- trans_rarefy$new(meco_amp_mean, alphadiv = "Pielou",  depth = c(0, 25, 100, 400, 800, 1600, 2400, 4800, 8000, 16000))
plot_amp_mean_rarefy_pielou <- meco_amp_mean_rarefy_pielou$plot_rarefy()

meco_amp_mean_rarefy <- meco_amp_mean_rarefy_sh$res_rarefy %>%
    left_join(meco_amp_mean_rarefy_simpson$res_rarefy, by = c("SampleID","seqnum")) %>% 
    left_join(meco_amp_mean_rarefy_pielou$res_rarefy, by = c("SampleID","seqnum")) %>%
    mutate(Category = "Mean")

################################
### END OF MEAN_AMP ANALYSIS ###
################################

## Create an Amp meco object for filtered ASV Counts
### Split the combined_amp to bact and arch tables
### Remove any Zeroes 
combined_zero_ASVs <- combined_amp %>%
  rowwise() %>%
  filter(!(sum(c_across(DMO.GOWN.C1.20200316:MDI.GOWN.E9.20200316), na.rm = TRUE) == 0)) %>%
  ungroup()
### MDI is Bact so we will remove any Archaea with counts giving back around 697 bact ASVs
MDI_combined_amp <- combined_amp %>%
    select(ASVs,Kingdom,Phylum,Class,Order,Family,Genus,starts_with("MDI")) %>%
    rename_with(
        ~ str_split(.x, "\\.", simplify = TRUE)[,3],
        .cols = -c(ASVs,Kingdom,Phylum,Class,Order,Family,Genus)
    ) %>%
    rowwise() %>%
    filter(!(sum(c_across(C1:E9), na.rm = TRUE) == 0) ) %>%
    ungroup() %>%
    filter(Kingdom == "Bacteria")

### DMO is Archaeal Primers so we will filter out any bacteria, leaving us 60 ASVs
DMO_combined_amp <- combined_amp %>%
    select(ASVs,Kingdom,Phylum,Class,Order,Family,Genus,starts_with("DMO")) %>%
    rename_with(
        ~ str_split(.x, "\\.", simplify = TRUE)[,3],
        .cols = -c(ASVs,Kingdom,Phylum,Class,Order,Family,Genus)
    )%>%
    rowwise() %>%
    filter(!(sum(c_across(C1:E9), na.rm = TRUE) == 0)) %>%
    ungroup() %>%
    filter(Kingdom == "Archaea")

### Here we will just combined the two and named it filtered with 757 ASVs
rearranged_filtered_combined_amp <- bind_rows(DMO_combined_amp, MDI_combined_amp) 

taxa_table_amp <- rearranged_filtered_combined_amp %>%
    select(ASVs,Kingdom,Phylum,Class,Order,Family,Genus) %>%
    column_to_rownames("ASVs")

otu_table_amp <- rearranged_filtered_combined_amp %>%
    select(-Kingdom,-Phylum,-Class,-Order,-Family,-Genus) %>%
    column_to_rownames("ASVs")

### Create the meco object for mean counts
meco_amp_filtered <- microtable$new(
    otu_table = otu_table_amp,
    tax_table = taxa_table_amp,
    sample_table = metadata
)
meco_amp_filtered$cal_alphadiv(measures=c("Observed","Shannon","Simpson","Pielou"))

### Identify total counts per sample
total_counts_amp <- colSums(meco_amp_filtered$otu_table)
min(total_counts_amp)
max(total_counts_amp)

meco_amp_filtered_rarefy_sh <- trans_rarefy$new(meco_amp_filtered, alphadiv = "Shannon",  depth = c(0, 25, 100, 400, 800, 1600, 2400, 4800, 8000, 16000))
plot_amp_filtered_rarefy_shannon <- meco_amp_filtered_rarefy_sh$plot_rarefy()

meco_amp_filtered_rarefy_simpson <- trans_rarefy$new(meco_amp_filtered, alphadiv = "Simpson",  depth = c(0, 25, 100, 400, 800, 1600, 2400, 4800, 8000, 16000))
plot_amp_filtered_rarefy_simpson <- meco_amp_filtered_rarefy_simpson$plot_rarefy()

meco_amp_filtered_rarefy_pielou <- trans_rarefy$new(meco_amp_filtered, alphadiv = "Pielou",  depth = c(0, 25, 100, 400, 800, 1600, 2400, 4800, 8000, 16000))
plot_amp_filtered_rarefy_pielou <- meco_amp_filtered_rarefy_pielou$plot_rarefy()

meco_amp_filtered_rarefy_os <- trans_rarefy$new(meco_amp_filtered, alphadiv = "Observed",  depth = c(0, 25, 100, 400, 800, 1600, 2400, 4800, 8000, 16000))
plot_amp_filtered_rarefy_os <- meco_amp_filtered_rarefy_os$plot_rarefy()

meco_amp_filtered_rarefy <- meco_amp_filtered_rarefy_sh$res_rarefy %>%
    left_join(meco_amp_filtered_rarefy_simpson$res_rarefy, by = c("SampleID","seqnum")) %>% 
    left_join(meco_amp_filtered_rarefy_pielou$res_rarefy, by = c("SampleID","seqnum")) %>% 
    left_join(meco_amp_filtered_rarefy_os$res_rarefy, by = c("SampleID","seqnum")) %>%
    mutate(Category = "Filtered")

####################################
### END OF Filtered_AMP ANALYSIS ###
####################################

# Processing the Metagenome Data and remove doubletons and singletons
undivided_combined_mtg <- combined_mtg %>%
    # select(ASVs,Domain,Phylum,Class,Order,Family,Genus,starts_with("MDI")) %>%
    rename_with(
        ~ str_split(.x, "\\.", simplify = TRUE)[,2],
        .cols = -c(ASVs,Domain,Phylum,Class,Order,Family,Genus)
    ) %>%
    rowwise() %>%
    filter(!(sum(c_across(C1:E9), na.rm = TRUE) <= 2) ) %>%
    ungroup()


### Making MECO_MTG object

taxa_table_mtg <- undivided_combined_mtg %>%
    select(ASVs,Domain,Phylum,Class,Order,Family,Genus) %>%
    column_to_rownames("ASVs")

otu_table_mtg <- undivided_combined_mtg %>%
    select(-Domain,-Phylum,-Class,-Order,-Family,-Genus) %>%
    column_to_rownames("ASVs")

### Create the meco object for mean counts
meco_mtg_filtered <- microtable$new(
    otu_table = otu_table_mtg,
    tax_table = taxa_table_mtg,
    sample_table = metadata
)
meco_mtg_filtered$cal_alphadiv(measures=c("Observed","Shannon","Simpson","Pielou"))

### Identify total counts per sample
total_counts_mtg <- colSums(meco_mtg_filtered$otu_table)
min(total_counts_mtg)
max(total_counts_mtg)

meco_mtg_filtered_rarefy_sh <- trans_rarefy$new(meco_mtg_filtered, alphadiv = "Shannon",  depth = c(0, 25, 100, 400, 800, 1600, 2400, 4800, 8000, 16000))
plot_mtg_filtered_rarefy_shannon <- meco_mtg_filtered_rarefy_sh$plot_rarefy()

meco_mtg_filtered_rarefy_simpson <- trans_rarefy$new(meco_mtg_filtered, alphadiv = "Simpson",  depth = c(0, 25, 100, 400, 800, 1600, 2400, 4800, 8000, 16000))
plot_mtg_filtered_rarefy_simpson <- meco_mtg_filtered_rarefy_simpson$plot_rarefy()

meco_mtg_filtered_rarefy_pielou <- trans_rarefy$new(meco_mtg_filtered, alphadiv = "Pielou",  depth = c(0, 25, 100, 400, 800, 1600, 2400, 4800, 8000, 16000))
plot_mtg_filtered_rarefy_pielou <- meco_mtg_filtered_rarefy_pielou$plot_rarefy()

meco_mtg_filtered_rarefy_os <- trans_rarefy$new(meco_mtg_filtered, alphadiv = "Observed",  depth = c(0, 25, 100, 400, 800, 1600, 2400, 4800, 8000, 16000))
plot_mtg_filtered_rarefy_os <- meco_mtg_filtered_rarefy_os$plot_rarefy()

meco_mtg_filtered_rarefy <- meco_mtg_filtered_rarefy_sh$res_rarefy %>%
    left_join(meco_mtg_filtered_rarefy_simpson$res_rarefy, by = c("SampleID","seqnum")) %>% 
    left_join(meco_mtg_filtered_rarefy_pielou$res_rarefy, by = c("SampleID","seqnum")) %>% 
    left_join(meco_mtg_filtered_rarefy_os$res_rarefy, by = c("SampleID","seqnum")) %>%
    mutate(Category = "Metagenome")

# Plot the different rarefaction curves
combined_rarefaction <- bind_rows(meco_amp_filtered_rarefy,meco_amp_mean_rarefy,meco_mtg_filtered_rarefy)

combined_rarefaction_long <- combined_rarefaction %>% 
  mutate(
    Category = factor(Category, levels = c("Filtered", "Mean", "Metagenome"))
  ) %>%
  pivot_longer(
    cols = c(Observed, Shannon, Pielou, Simpson),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = c("Observed", "Shannon", "Pielou", "Simpson"))
  )

plot_combined_rarefaction <- ggplot(combined_rarefaction_long, 
                                aes(x = seqnum, y = Value, group = SampleID)
                                ) +
                                geom_line(aes(color = SampleID, linetype = Category),
            linewidth = 0.7, alpha = 0.9) +
  facet_grid(Metric ~ Category, scales = "free_y") +
  scale_linetype_manual(values = c(
    Filtered   = "solid",
    Mean       = "dashed",
    Metagenome = "dotted"
  )) +
  labs(
    x = "Sequences subsampled",
    y = "Diversity / Richness",
    color = "Sample",
    linetype = "Category",
    title = "Rarefaction curves by category and diversity metric"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

plot_combined_rarefaction
ggsave(paste0(outpath,"plot_combined_rarefaction.png"), plot_combined_rarefaction)

######################################## END OF ALPHA DIVERSITY#############################################

######################################
# BETA DIVERSITY PAIRWISE COMPARISON #
######################################

meco_mtg_filtered$cal_betadiv(method="bray")
meco_amp_filtered$cal_betadiv(method="bray")
meco_amp_mean$cal_betadiv(method="bray")

beta_mtg <- trans_beta$new(meco_mtg_filtered, measure = "bray")
beta_amp_filtered <- trans_beta$new(meco_amp_filtered, measure = "bray")
beta_amp_mean <- trans_beta$new(meco_amp_mean, measure = "bray")

beta_mtg_long <- beta_mtg$use_matrix %>%
     as.data.frame() %>%
  rownames_to_column(var = "Row") %>%
  pivot_longer(
    cols = -Row,
    names_to = "Column",
    values_to = "Value"
  ) %>%
  mutate(Category = "Metagenome")

beta_amp_filtered_long <- beta_amp_filtered$use_matrix %>%
     as.data.frame() %>%
  rownames_to_column(var = "Row") %>%
  pivot_longer(
    cols = -Row,
    names_to = "Column",
    values_to = "Value"
  ) %>%
  mutate(Category = "Filtered")


beta_amp_mean_long <- beta_amp_mean$use_matrix %>%
     as.data.frame() %>%
  rownames_to_column(var = "Row") %>%
  pivot_longer(
    cols = -Row,
    names_to = "Column",
    values_to = "Value"
  ) %>%
  mutate(Category = "Mean")

beta_pairwise_distances <- bind_rows(beta_amp_mean_long, beta_amp_filtered_long, beta_mtg_long)

wide_df <- beta_pairwise_distances %>%
  pivot_wider(
    names_from = Category,
    values_from = Value
  ) %>%
  # remove self-pairs if present
  filter(Row != Column)

beta_mean_filtered <- ggplot(wide_df, aes(x = Mean, y = Filtered)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, color = "black") +
  labs(
    x = "Mean distances",
    y = "Filtered distances",
    title = "Pairwise Beta Diversity: Mean vs Filtered"
  ) +
  theme_bw(base_size = 12)

beta_mean_metagenome <- ggplot(wide_df, aes(x = Mean, y = Metagenome)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, color = "black") +
  labs(
    x = "Mean distances",
    y = "Metagenome distances",
    title = "Pairwise Beta Diversity: Mean vs Metagenome"
  ) +
  theme_bw(base_size = 12)

beta_filtered_metagenome <- ggplot(wide_df, aes(x = Filtered, y = Metagenome)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, color = "black") +
  labs(
    x = "Filtered distances",
    y = "Metagenome distances",
    title = "Pairwise Beta Diversity: Filtered vs Metagenome"
  ) +
  theme_bw(base_size = 12)


library(patchwork)

beta_pairwise <-  beta_filtered_metagenome + beta_mean_filtered + beta_mean_metagenome
ggsave("gownproject/P2-Strain_Heterogeneity/100_AmpliconMetagenomeComparison/10_plot_beta_pairwise.png")
