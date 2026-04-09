library(tidyverse)

# Input data
singlem_profile_path <- "gownproject/P3-Mesocosm/100_denitrificationpotential_MESO/30_taxonomicabundance/03_MESO25/03_MESO25_combined-profiles.summarise.wextras.tsv"
output_path <- "gownproject/P3-Mesocosm/100_denitrificationpotential_MESO/90_Figures"
denitrificaion_path <- "/home/gam/github/gownproject/P3-Mesocosm/900_wetlab/Results/DenitrificationRate/DenitrificationRate_SED_Individual_Slopes.csv"

dnr_df <- read_csv(denitrificaion_path)

df <-read_tsv(singlem_profile_path) %>%
  select(sample, taxonomy, level, relative_abundance) %>%
  filter(!is.na(relative_abundance)) -> df


#----------------------------
# 1) Split taxonomy into ranks
#----------------------------
df2 <- df %>%
  mutate(
    taxonomy = str_remove(taxonomy, "^Root;\\s*")
  ) %>%
  separate(
    taxonomy,
    into = c("domain","phylum","class","order","family","genus","species"),
    sep = ";\\s*",
    fill = "right"
  )

#----------------------------
# 2) Create parent keys
#----------------------------
df2 <- df2 %>%
  mutate(
    parent_key = case_when(
      level == "species" ~ paste(sample, domain, phylum, class, order, family, genus, sep="|"),
      level == "genus"   ~ paste(sample, domain, phylum, class, order, family, sep="|"),
      level == "family"  ~ paste(sample, domain, phylum, class, order, sep="|"),
      level == "order"   ~ paste(sample, domain, phylum, class, sep="|"),
      level == "class"   ~ paste(sample, domain, phylum, sep="|"),
      level == "phylum"  ~ paste(sample, domain, sep="|"),
      level == "domain"  ~ paste(sample, sep="|"),
      TRUE ~ NA_character_
    ),
    self_key = case_when(
      level == "species" ~ paste(sample, domain, phylum, class, order, family, genus, species, sep="|"),
      level == "genus"   ~ paste(sample, domain, phylum, class, order, family, genus, sep="|"),
      level == "family"  ~ paste(sample, domain, phylum, class, order, family, sep="|"),
      level == "order"   ~ paste(sample, domain, phylum, class, order, sep="|"),
      level == "class"   ~ paste(sample, domain, phylum, class, sep="|"),
      level == "phylum"  ~ paste(sample, domain, phylum, sep="|"),
      level == "domain"  ~ paste(sample, domain, sep="|"),
      TRUE ~ paste(sample, sep="|")
    )
  )

#----------------------------
# 3) Sum children per parent
#----------------------------
child_sum <- df2 %>%
  group_by(parent_key) %>%
  summarise(child_abundance = sum(relative_abundance), .groups="drop")

#----------------------------
# 4) Join and compute exclusive abundance
#----------------------------
df_fixed <- df2 %>%
  left_join(child_sum, by = c("self_key" = "parent_key")) %>%
  mutate(
    child_abundance = replace_na(child_abundance, 0),
    exclusive_abundance = relative_abundance - child_abundance,
    exclusive_abundance = pmax(exclusive_abundance, 0) # avoid negatives
  )


# ----------------------------
# 5) Finalize output
# ----------------------------
tax_cols <-c("domain","phylum","class","order","family","genus","species")

final_df <- df_fixed %>%
  select(sample, all_of(tax_cols), level, exclusive_abundance) %>%
  rename(relative_abundance = exclusive_abundance) %>%
  filter(level != "root") %>%
  mutate(TAX_ID = paste0("TAX_", row_number()))

singlem_profile <- final_df %>%
  pivot_wider(
    id_cols = c(all_of(tax_cols), level, TAX_ID),
    names_from = sample,
    values_from = relative_abundance,
    values_fill = 0
  )

# Fixed taxonomic names
## Functions
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

tax_filled <- singlem_profile %>%
    mutate(domain = ifelse(domain == "Unassigned", "d__Unassigned", domain)) %>%
    mutate(
        across(
        domain:species,
        ~ {
            x <- replace_na(.x, "")
            core <- str_remove(x, "^[a-z]__")
            if_else(core %in% c("", "unclassified"), "unclassified", x)
        }
        )
    )

tax_filled2 <- tax_filled
tax_filled2[tax_cols] <- t(apply(tax_filled2[tax_cols],1, fill_from_parent))

singlem_df <- tax_filled2 %>%
  rename_with(~ str_extract(.x, "(CO|SA|SS|SH)\\d+"), 
              .cols = -c(level, TAX_ID, all_of(tax_cols)))

singlem_df_filtered <- singlem_df %>%
  filter(rowSums(across(CO1:SS3)) > 0)

# ----------------------------
# 6) Create Microeco object
# ----------------------------
library(microeco)
tax_table <- singlem_df_filtered %>%
  select(all_of(tax_cols), TAX_ID) %>%
  column_to_rownames("TAX_ID")

otu_table <- singlem_df_filtered %>%
  select(-all_of(tax_cols), -level) %>%
  column_to_rownames("TAX_ID")

 metadata <- levels(factor(final_df$sample)) %>%
  as.data.frame() %>%
  mutate(SampleID = str_extract(., "(CO|SA|SS|SH)\\d+"),
         Day = str_extract(., "D\\d+"),
         Setup = case_when(
           str_starts(SampleID, "CO") ~ "Coal",
           str_starts(SampleID, "SA") ~ "Sand",
           str_starts(SampleID, "SS") ~ "Sandstone",
           str_starts(SampleID, "SH") ~ "Shale",
           TRUE ~ NA_character_
         ),
         Replicate = str_extract(SampleID, "\\d+$")
  ) %>%
  left_join(
    dnr_df %>% 
      select(Setup, Replicate, slope) %>% 
      mutate(Replicate = as.character(Replicate)), 
    by = c("Setup", "Replicate")
  ) %>%
  select(Day, Setup, Replicate, slope, SampleID) %>%
  rename(Denitrification_Slope = slope) %>%
  column_to_rownames("SampleID")

# ----------------------------
# Microeco MAIN
# ----------------------------
meso_meco <- microtable$new(
  otu_table = otu_table,
  tax_table = tax_table,
  sample_table = metadata
)

meso_meco$cal_abund()
#cols <- c("#6ea577", "#D7B18C", "#F0A390", "#AAB0D4", "#D0A5D8", "#ECA8BE", "#BDD676", "#D7E151", "#F28E2B", "#F6D383", "#DFC9AB", "#C2C2C2", "grey70")
cols <- c(
  # original pastel (12)
  "#6ea577", "#D7B18C", "#F0A390", "#AAB0D4", "#D0A5D8", "#ECA8BE",
  "#BDD676", "#D7E151", "#F28E2B", "#F6D383", "#DFC9AB", "#C2C2C2",
  
  # new (more distinct / less pastel)
  "#2E7D32",  # deeper green
  "#C97C30",  # burnt orange
  "#D55E5E",  # muted red
  "#5B6FB8",  # stronger blue
  "#9A4FBF",  # deeper purple
  "#C94A7C"   # stronger pink
)
meso_abund <- trans_abund$new(meso_meco, taxrank = "genus", ntaxa = 18)
plot_abund_sed <- meso_abund$plot_bar(
  facet = "Setup",
  xtext_keep = FALSE,
  legend_text_italic = TRUE,
  others_color = "grey40",
  color_values = cols) 

plot_abund_sed

# -----------------------------
# Environmental correlation
# -----------------------------
meso_env <- trans_env$new(meso_meco, env_cols = c(2,4))
meso_env$cal_diff(group = "Setup", method = "wilcox")
meso_env$res_diff

meso_env$cal_ordination(method = "RDA", dist_method = "bray", taxa_level="genus")
meso_env$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE,  max_perc_env = 1.5, max_perc_tax = 1.0)
plot_env <- meso_env$plot_ordination(plot_color = "Setup")
plot_env

# ggsave(
#   plot = plot_abund_sed,
#   filename = file.path(output_path, "Singlem_Abundance_Barplot_Genus.png"),
#   width = 10, height = 6, dpi = 300
# )

# ggsave(
#   plot = plot_env,
#   filename = file.path(output_path, "Singlem_RDA_Ordination_Genus.png"),
#   width = 8, height = 6, dpi = 300
# )
