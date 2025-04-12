# comparison_fqdept.R
library(tidyverse)
library(microeco)
library(ggpubr)
library(ggsci)

# import all tsv species file
input <- '/home/glbcabria/Workbench/P1/subset_fastq/singlem'
files <- list.files(input, pattern = "_relabund-species\\.tsv", full.names = TRUE)

for (file in files) {
    filename <- basename(file)
    base_name <- str_split_i(filename,pattern = "_relabund",1)
    var_name <- make.names(base_name)
    
    replicate <- str_split_i(var_name, patter = "\\.rep", -1)
    depth <- str_split_i(var_name, patter = "\\.", -2)
    
    df <- read_tsv(file) %>%
      rename_with(~ "Fraction", .cols =2) %>%
      mutate(Depth = depth) %>%
      mutate(Replicate = replicate)
    
    assign(var_name, df, envir = .GlobalEnv)
}

df_names <- ls(pattern ="rep\\d+$")
df_list <- lapply(df_names, get)
combined_df <- bind_rows(df_list, .id = "source")

relabund_matrix <- combined_df %>%
  unite(col = 'Col_ID', Depth:Replicate, sep = '_',) %>%
  filter(!grepl('unsampled', Col_ID)) %>%
  select(taxonomy, Fraction, Col_ID) %>%
  pivot_wider(names_from = Col_ID,
              values_from = Fraction) %>%
  mutate(OTU_ID = paste0("OTU_", row_number() ) ) %>%
  column_to_rownames("OTU_ID")

otu_table <-relabund_matrix %>%
  select(-taxonomy) %>%
  mutate(across (
    everything(),
    ~ifelse(is.na(.), 0, .)
  )
  )

taxa_table <- relabund_matrix %>%
  select(taxonomy) %>%
  separate(taxonomy, into = c('Root','Kingdom','Phylum','Class','Order','Family','Genus','Species'), sep = ';') %>%
  tidy_taxonomy() %>%
  mutate( across(
    everything(), 
    ~ ifelse(Root == 'r__', gsub('__','__unassigned',.), .)
      )
    ) 

sample_table <- data.frame(Cols = colnames(relabund_matrix)[-1]) %>%
  separate(Cols, into = c('Subsample_Size','Replicate'), remove = FALSE) %>%
  column_to_rownames('Cols')

#####################
# Microeco Analysis #
#####################

meco <- microtable$new(
  otu_table = otu_table,
  tax_table = taxa_table,
  sample_table = sample_table
)

meco$sample_table$Subsample_Size <- factor(
  meco$sample_table$Subsample_Size,
  levels = c('unsampled','20M','10M','5M','1M','500K')
)

meco$cal_abund()
meco$cal_alphadiv()
# t1 <-trans_venn$new(dataset = meco)
# g1 <- t1$plot_bar(left_plot = TRUE, bottom_height = 0.5, left_width = 0.15, up_bar_fill = "grey50", left_bar_fill = "grey50", bottom_point_color = "black")

t2 <- trans_alpha$new(dataset = meco, group = "Subsample_Size")
t2$cal_diff(method = "wilcox")
Shannon <- t2$plot_alpha(measure = "Shannon")
Shannon$layers[[2]]<-NULL
Shannon + geom_jitter(width = 0.2, size =2, alpha = 0.6) + scale_color_jco() +
  stat_compare_means(method = "wilcox", comparisons = list(c('20M','1M')), label = "p.signif")


Simpson <- t2$plot_alpha(plot_type = "ggboxplot", measure = "Simpson") #, y_increase = 0.2, add_line = TRUE, line_type = 2, line_alpha = 0.5, errorbar_width = 0.1 )
InvSimpson <- t2$plot_alpha(plot_type = "ggboxplot", measure = "InvSimpson") #, y_increase = 0.2, add_line = TRUE, line_type = 2, line_alpha = 0.5, errorbar_width = 0.1 )
Pielou <-t2$plot_alpha(plot_type = "ggboxplot", measure = "Pielou") #, y_increase = 0.2, add_line = TRUE, line_type = 2, line_alpha = 0.5, errorbar_width = 0.1 )
Coverage <- t2$plot_alpha(plot_type = "errorbar", measure = "Coverage") #, y_increase = 0.2, add_line = TRUE, line_type = 2, line_alpha = 0.5, errorbar_width = 0.1 )

# adjust_label_y <- function(p, y, z=0.1) {
#   p + stat_compare_means(label = "p.signif", label.y = y, step.increase = z)
# }
# 
# Shannon <- adjust_label_y(Shannon, y = 6, z= 0.001)

Shannon

ggarrange(Shannon, Simpson, InvSimpson, Pielou, Coverage,
          ncol = 3, nrow = 2,
          labels = c("A", "B", "C", "D", "E")
)
