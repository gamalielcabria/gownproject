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
    as.data.frame() %>%
    mutate(pH = c(8.5, 8.5, 8.5, 10.2, 10.2,10.2, 8.5, 8.5),
           Temperature = c(21,21,21,21,21,21,30,30),
           Alkalinity = c(0.24,0.24,0.24,0.5,0.5,0.5,0.24,0.24)
    )



meco <- microtable$new(otu_table = otu_table, sample_table = sample_metadata)
meco$cal_betadiv(measure = "bray")

beta_diversity <- trans_beta$new(dataset = meco, group = "Site", measure = "bray")
beta_diversity$cal_ordination(method = "NMDS")

plot_bd <- beta_diversity$plot_ordination(plot_color = "Site", plot_shape = "Timepoint")
plot_bd

library(ggokabeito)

nmds_beautiful <- plot_bd +
    geom_point(size = 8) +
    scale_color_okabe_ito() +
    theme_minimal() +
    theme(
        aspect.ratio = 1,
        text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
    ) +
    labs(
        title = "NMDS of Microbial Communities",
        x = "NMDS1",
        y = "NMDS2",
        color = "Site",
        shape = "Timepoint"
    ) +
    guides(
        color = guide_legend(override.aes = list(size=5)),
        shape = guide_legend(override.aes = list(size=5))
    )

nmds_beautiful

ggsave(
    nmds_beautiful,
    filename = "gownproject/Script/utility_scripts/all_nmds.png",
    width = 8,
    height = 6,
    dpi = 300
)

bray_curtis <- beta_diversity$use_matrix
write.csv(bray_curtis, file = "gownproject/Script/utility_scripts/bray_curtis.csv")

# Environmental Fitting
env_fit <- trans_env$new(dataset = meco, env_cols = 3:5)
env_fit$cal_ordination(method = "CCA", use_measure = "bray")
env_fit$trans_ordination()
plot_ef <- env_fit$plot_ordination(plot_color = "Site", plot_shape = "Timepoint")
plot_ef
