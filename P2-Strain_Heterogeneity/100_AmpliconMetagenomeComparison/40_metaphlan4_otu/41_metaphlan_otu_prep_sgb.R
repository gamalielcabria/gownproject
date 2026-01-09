# Script for preparing the qiime taxonomy and otu table inputs
library(tidyverse)
library(microeco)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

# Set Working Directory~
#setwd("/home/glbcabria/Workbench")

# Import OTU Counts, taxa, and metadata
## File Path of inputs
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
#     separate(otu_id, into = tax_cols, sep = "\\|", fill = "right", extra = "drop") %>%
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
    separate(otu_id, into = tax_cols[1:7], sep = "\\|", fill = "right", extra = "drop") %>%
    rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
        .cols = -c(Kingdom,Phylum,Class,Order,Family,Genus,Species)
    ) %>%
    select(
        all_of(tax_cols[1:7]),
        sort(setdiff(colnames(.),tax_cols))
    ) %>%
    mutate(OTUs = ifelse(grepl("k__Bacteria", Kingdom), paste0("bOTU_", row_number()), paste0("aOTU_", row_number())))

sgb_genus_otu_tax <- assigned_otu_tax_sgb %>%
    filter(grepl("g__", otu_id) & !grepl("s__", otu_id)) %>%
    separate(otu_id, into = tax_cols[1:6], sep = "\\|", fill = "right", extra = "drop") %>%
    rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
        .cols = -c(Kingdom,Phylum,Class,Order,Family,Genus)
    ) %>%
    select(
        all_of(tax_cols[1:6]),
        sort(setdiff(colnames(.),tax_cols))
    ) %>%
    mutate(OTUs = ifelse(grepl("k__Bacteria", Kingdom), paste0("bOTU_", row_number()), paste0("aOTU_", row_number())))

sgb_family_otu_tax <- assigned_otu_tax_sgb %>%
    filter(grepl("f__", otu_id) & !grepl("g__", otu_id) & !grepl("s__", otu_id)) %>%
    separate(otu_id, into = tax_cols[1:5], sep = "\\|", fill = "right", extra = "drop") %>%
    rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
        .cols = -c(Kingdom,Phylum,Class,Order,Family)
    ) %>%
    select(
        all_of(tax_cols[1:5]),
        sort(setdiff(colnames(.),tax_cols))
    ) %>%
    mutate(OTUs = ifelse(grepl("k__Bacteria", Kingdom), paste0("bOTU_", row_number()), paste0("aOTU_", row_number())))

sgb_phylum_otu_tax <- assigned_otu_tax_sgb %>%
    filter(grepl("p__", otu_id) & !grepl("c__", otu_id)) %>%
    separate(otu_id, into = tax_cols[1:2], sep = "\\|", fill = "right", extra = "drop") %>%
    rename_with(~ str_split(.x, "-", simplify = TRUE)[,2],
        .cols = -c(Kingdom,Phylum)
    ) %>%
    select(
        all_of(tax_cols[1:2]),
        sort(setdiff(colnames(.),tax_cols))
    ) %>%
    mutate(OTUs = ifelse(grepl("k__Bacteria", Kingdom), paste0("bOTU_", row_number()), paste0("aOTU_", row_number())))

## Save as csv files
# write.csv(sgb_species_otu_tax, file = paste0(outpath, "/mp4-species_bact_arch_combined_otu.csv"))
# write.csv(sgb_genus_otu_tax, file = paste0(outpath, "/mp4-genus_bact_arch_combined_otu.csv"))
# write.csv(sgb_family_otu_tax, file = paste0(outpath, "/mp4-family_bact_arch_combined_otu.csv"))
# write.csv(sgb_phylum_otu_tax, file = paste0(outpath, "/mp4-phylum_bact_arch_combined_otu.csv"))

#################
# END OF SCRIPT #
#################

# ## Function to create rank tables with unassigned counts
# make_rank_tables_with_unassigned <- function(df, id_col = "otu_id") {
#   df <- as_tibble(df)
#   stopifnot(id_col %in% names(df))

#   sample_cols <- setdiff(names(df), id_col)

#   ranks <- c(
#     kingdom = "k__", phylum = "p__", class = "c__", order = "o__",
#     family  = "f__", genus  = "g__", species = "s__"
#   )

#   get_rank_vec <- function(x, prefix) {
#     # last matching token in path (or NA)
#     hits <- str_extract_all(x, "\\b[a-z]__[^|]+") # all rank tokens in order
#     map_chr(hits, \(v) {
#       v2 <- v[startsWith(v, prefix)]
#       if (length(v2) == 0) NA_character_ else v2[length(v2)]
#     })
#   }

#   x <- df[[id_col]]
#   rank_cols <- imap_dfc(ranks, \(pref, nm) tibble(!!nm := get_rank_vec(x, pref)))
#   dat <- bind_cols(df, rank_cols)

#   # ---- find leaves (avoid double-counting rollups already present) ----
#   parent_paths <- dat %>%
#     transmute(parent = if_else(str_detect(.data[[id_col]], "\\|"),
#                               str_replace(.data[[id_col]], "\\|[^|]+$", ""),
#                               NA_character_)) %>%
#     filter(!is.na(parent)) %>%
#     distinct(parent) %>%
#     pull(parent)

#   leaves <- dat %>% filter(!( .data[[id_col]] %in% parent_paths ))

#   # ---- aggregate + unassigned helper (child within parent) ----
#   agg_rank <- function(child_rank, parent_rank = NULL) {
#     child_prefix <- ranks[[child_rank]]

#     # counts at the child rank from leaves
#     child_tbl <- leaves %>%
#       filter(!is.na(.data[[child_rank]])) %>%
#       group_by(rank_id = .data[[child_rank]]) %>%
#       summarise(across(all_of(sample_cols), sum), .groups = "drop")

#     # root (kingdom): also include anything with NA kingdom
#     if (is.null(parent_rank)) {
#       unassigned_root <- leaves %>%
#         filter(is.na(.data[[child_rank]])) %>%
#         summarise(across(all_of(sample_cols), sum)) %>%
#         mutate(rank_id = paste0(child_prefix, "UNASSIGNED")) %>%
#         select(rank_id, all_of(sample_cols))

#       return(bind_rows(child_tbl, unassigned_root) %>% arrange(rank_id))
#     }

#     # totals at parent from leaves
#     parent_totals <- leaves %>%
#       filter(!is.na(.data[[parent_rank]])) %>%
#       group_by(parent_id = .data[[parent_rank]]) %>%
#       summarise(across(all_of(sample_cols), sum), .groups = "drop")

#     # sum of leaves that reach the child rank within each parent
#     reached_child <- leaves %>%
#       filter(!is.na(.data[[parent_rank]]), !is.na(.data[[child_rank]])) %>%
#       group_by(parent_id = .data[[parent_rank]]) %>%
#       summarise(across(all_of(sample_cols), sum), .groups = "drop")

#     joined <- parent_totals %>%
#       left_join(reached_child, by = "parent_id", suffix = c("_parent", "_reached"))

#     # matrix subtraction for unassigned = parent - reached
#     parent_mat  <- as.matrix(joined %>% select(ends_with("_parent")))
#     reached_mat <- as.matrix(joined %>% select(ends_with("_reached"))) %>% { replace(., is.na(.), 0) }
#     diff_mat <- parent_mat - reached_mat
#     diff_mat[diff_mat < 0] <- 0  # guard tiny negatives from rounding

#     # name columns back to sample names
#     colnames(diff_mat) <- sample_cols

#     unassigned_child <- bind_cols(
#       tibble(rank_id = paste0(child_prefix, "UNASSIGNED_in_", joined$parent_id)),
#       as_tibble(diff_mat)
#     )

#     bind_rows(child_tbl, unassigned_child) %>% arrange(rank_id)
#   }

#   list(
#     kingdom = agg_rank("kingdom", NULL),
#     phylum  = agg_rank("phylum",  "kingdom"),
#     class   = agg_rank("class",   "phylum"),
#     order   = agg_rank("order",   "class"),
#     family  = agg_rank("family",  "order"),
#     genus   = agg_rank("genus",   "family"),
#     species = agg_rank("species", "genus")
#   )
# }

# # ---- run ----
# rank_tables <- make_rank_tables_with_unassigned(otu_tax_sgb_mdf, id_col = "otu_id")

# # examples:
# rank_tables$species  # includes s__... plus s__UNASSIGNED_in_g__...
# rank_tables$phylum   # includes p__... plus p__UNASSIGNED_in_k__...
