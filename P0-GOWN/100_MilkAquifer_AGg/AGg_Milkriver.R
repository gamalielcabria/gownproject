library("tidyverse")
library("microeco")
library("janitor")

# read in data 
setwd("gownproject/P0-GOWN/100_MilkAquifer_AGg/")
tax_raw <- read_csv("taxa.csv", col_names = TRUE) %>%
    setNames(c("Sequence","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
    
otu_raw <- read_csv("otu_table.csv", col_names = FALSE) %>%
    t() %>%
    as.data.frame() %>%
    row_to_names(row_number =1 )  %>%
    setNames(c("Sequence","GOWN211","GOWN213","GOWN212")) %>%
    left_join(tax_raw, by = "Sequence") %>%
    mutate(ASVs = paste0("ASV", 1:nrow(.))) %>%
    #rownames_to_column("old") %>%
    column_to_rownames("ASVs") %>%
    select(-Sequence)

# create tax table
tax <- otu_raw %>%
    select(Kingdom:Species) %>%
    tidy_taxonomy()  %>%
  mutate(across(everything(), ~ ifelse(grepl("__$", .x), NA, .x)))

tax_filled <- t(apply(tax, 1, function(row) {
  last <- NA
  for (i in seq_along(row)) {
    if (!is.na(row[i])) last <- row[i]
    else row[i] <- last
  }
  row
})) %>% as.data.frame()

otu <- otu_raw %>%
    select(GOWN211:GOWN212) %>%
    as.data.frame() %>%
    mutate(across(everything(), as.numeric))

meco_milkriver <- microtable$new(otu = otu, tax = tax_filled)

abund_milkriver <- trans_abund$new(meco_milkriver, taxrank = "Genus",ntaxa=20)
plot_abund <- abund_milkriver$plot_bar() +
    ggtitle("Milk River Aquifer Abundnace - Top 20 Genera") +
    ylab("Relative Abundance (%)") +
    xlab("GOWN Wells") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 20, hjust = 0.5),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)) +
  guides(
    color = guide_legend(ncol = 1),
    fill  = guide_legend(ncol = 1),
    shape = guide_legend(ncol = 1),
    linetype = guide_legend(ncol = 1)
  ) 

ggsave("AGg_MilkRiver_Abundance_Top20Genus.png", plot = plot_abund, width = 10, height = 8, dpi = 300)
