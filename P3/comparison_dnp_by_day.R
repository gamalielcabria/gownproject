# comparison_dnp_by_day.R
library("tidyverse")
library(scales)
library(ggpubr)
library(rstatix)
library(patchwork)
library(stringr)
library(ggsci)

input_d2 <- '/home/glbcabria/Workbench/P3/Results/GCPhotos/D2_stat_gc_norm.csv'
input_d43 <- '/home/glbcabria/Workbench/P3/Results/GCPhotos/D43_stat_gc_norm.csv'

stat_gc_norm_d2 <- read_csv(input_d2, col_names = TRUE) %>%
  mutate(DAY = 2)
stat_gc_norm_d43 <- read_csv(input_d43, col_names = TRUE)%>%
  mutate(DAY = 43)

comparison_gc_norm <- rbind(stat_gc_norm_d2, stat_gc_norm_d43)

stat_plot_norm <- ggboxplot(comparison_gc_norm,
                            x = 'DAY',
                            y = 'Normalized.Rate',
                            color = 'Rock',
                            add = c("jitter", "mean_sd"),
                            palette = "jco") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +    
  stat_compare_means(comparisons = list(c("2","43")),
                     method = "wilcox.test",
                     paired = FALSE, #hide.ns = TRUE,
                     label="p.signif")+
  #ylab("Denitrification Rate\n(nmol N2O per hour per cm^2)" ) +
  ylab( expression("Denitrification Rate"~"\n(" * nmol~N[2]*O~"/hour/" * cm^2 * ")") ) +
    facet_grid(~Rock) +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.placement = "outside",
    strip.background = element_rect(fill="grey90"),
    strip.text = element_text(face = "bold", size = 18),
    #strip.switch.pad.wrap = unit(0.2, "cm"),
    legend.position = "none"
  )

stat_plot_norm


# RUN THIS AFTER RUNNING GC_data_visualisation.R
combined_plot<- plot_DNP_agg_norm_rep + stat_plot_norm + stat_plot2 + 
  plot_annotation(tag_levels = "A", ) & 
  theme(plot.tag = element_text(size = 20) )
combined_plot
ggsave(plot = combined_plot, filename = "~/Workbench/P3/Results/combined_plot.png",
       units = c('px'),
       width = 3000, height = 1000, dpi = 120)
