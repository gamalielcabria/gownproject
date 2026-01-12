library(tidyverse)
library(microeco)

getwd()

input<-"gownproject/Script/utility_scripts/miseq_all.shared.txt"

mothur <- read_tsv(input, col_names = TRUE)

otu_table <- mothur %>%
    filter(label == "0.03") %>%
    select(-label, -numOtus) %>%
    column_to_rownames(var = "Group") %>%
    t() %>%
    as.data.frame()

sample_metadata <- colnames(otu_table) %>%
    as.data.frame() %>%
    rename(ID = ".") %>%
    separate( ID, into = c("Site", "Timepoint"), remove = FALSE) %>%
    column_to_rownames(var = "ID") %>%
    as.data.frame()

meco <- microtable$new(otu_table = otu_table, sample_table = sample_metadata)
meco$cal_betadiv(measure = "bray")

beta_diversity <- trans_beta$new(dataset = meco, group = "Site", measure = "bray")
beta_diversity$cal_ordination(method = "PCoA")

plot_bd <- beta_diversity$plot_ordination(plot_color = "Site", plot_shape = "Timepoint")
plot_bd
