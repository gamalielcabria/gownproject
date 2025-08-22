# SingleM Summarise to Plot
# This is a script to analyse the results of SingleM summarise
# with options '--output-taxonomic-profile-with-extras'

########################### Libraries#############################
library(tidyverse)
library(microeco)

############################ Inputs ##############################
# Input files has file extension: '.summarise.wextras.tsv' 
# but the extension can change can change
# Modify inputs based on your preference
inputfolder <- '/home/glbcabria/Workbench/P2/singlem/summary'
output <- 'singlem_GOWN2022'
file_ext <- "\\.summarise\\.wextras\\.tsv$" #Change if different extension
setwd(inputfolder)

############################# Main ###############################
# File Location
files <- list.files(
            pattern=file_ext, 
            ignore.case = TRUE
            )

# Create empty dataframe to store results
otus_df <- data.frame(
            sample = character(),
            coverage = numeric(),
            full_coverage = numeric(),
            relabund = numeric(),
            level = character(),
            taxonomy = character(),
            stringsAsFactors = FALSE
            )

# Read through files
for (file in files) {
    df <- read.csv(file, header = TRUE, sep = "\t")
    otus_df <- rbind(otus_df, df)
}

otu_taxa_table <- otus_df %>%
    mutate(sample = str_split_i(sample, '\\.', 1)) %>%
    select(sample, coverage, taxonomy) %>%
    pivot_wider(
        names_from = sample,
        values_from = coverage
    ) %>%
    mutate(taxonomy = gsub('Root; ','',taxonomy)) %>%
    mutate(taxonomy = gsub('Root','d__unassigned',taxonomy)) %>%
    mutate(id = paste0("OTU_",row_number()))

taxa_table <- otu_taxa_table %>%
    select(id, taxonomy) %>%
    column_to_rownames(var = 'id') %>%
    separate(
        col = taxonomy,  # replace with your column name
        into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus","Species"),
        sep = "; ",      # note the space after ;
        fill = "right",  # if some rows have fewer levels
        remove = FALSE   # keep original column, set to TRUE to drop it
    ) %>%
    tidy_taxonomy()

otu_table <- otu_taxa_table %>%
    select(-taxonomy) %>%
    column_to_rownames(var = 'id') %>%
    mutate(across(everything(), ~replace_na(.,0)))

######################### Write Output ############################
write.csv(file = paste0('otu_table_',output,'.csv'), otu_table)
write.csv(file = paste0('taxa_table_',output,'.csv'), otu_table)
