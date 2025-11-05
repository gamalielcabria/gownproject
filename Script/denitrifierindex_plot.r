# Denitrifier Abundance vs Index
library(tidyverse)
library(ggplot2)
library(readxl)
library(scales)
library(patchwork)

# Set WD and load inpuit
workdir <- '/home/glbcabria/Workbench/P0/'
setwd(workdir)
input <- 'DenitrificationPotential_1.xlsx'
input2 <- 'denitrification_list.txt'

dp_df <- read_excel(input, 
            sheet = 'Sheet1', 
            col_names = TRUE)

dp_metadata <- read_delim(input2,
                delim = '\t',
                col_names = TRUE
                )

headers <- c('SampleID',
            'Denitrifier_Total_Ratio',
            'Denitrifier_Rel_Abund',
            'Index', # %D/%O
            'Sample',
            'Site',
            'Nitrate' # mg/L
            )
            
colnames(dp_df) <- headers

dp_df2 <- dp_df %>%
            filter(Nitrate != 'NA') %>%
            mutate(Nitrate = as.numeric(Nitrate),
                    Index = ifelse(Index == '?', 0, Index),
                    Index = as.numeric(Index)
            ) %>%
            left_join(dp_metadata, by=c('Site'='Well')) %>%
            select(
                c(Denitrifier_Rel_Abund,Index,Nitrate,Site,
                `GOWN Lithology`, `Well Depth`, `Aquifer Type`
                )
            ) %>%
            rename(
                GOWN_Lithology = `GOWN Lithology`,
                Well_Depth = `Well Depth`,
                Aquifer_Type = `Aquifer Type`
            ) %>%
            mutate( Well_Depth = as.numeric(str_extract(Well_Depth, "\\d+\\.?\\d*")) ) %>%
            mutate(
                Shale     = str_detect(GOWN_Lithology, regex("Shale", ignore_case = TRUE)),
                Sand      = str_detect(GOWN_Lithology, regex("Sand\\b", ignore_case = TRUE)),
                Sandstone = str_detect(GOWN_Lithology, regex("Sandstone", ignore_case = TRUE)),
                Coal      = str_detect(GOWN_Lithology, regex("Coal", ignore_case = TRUE)),
                Gravel     = str_detect(GOWN_Lithology, regex("Gravel", ignore_case = TRUE)),
                Clay      = str_detect(GOWN_Lithology, regex("Clay", ignore_case = TRUE)),
                Limestone = str_detect(GOWN_Lithology, regex("Limestone", ignore_case = TRUE)),
                Unknown = str_detect(GOWN_Lithology, regex("unknown", ignore_case = TRUE))
                ) %>%
            mutate(Site = str_remove(Site, "\\s*\\([^\\)]+\\)"))

dp_df3 <- dp_df2 %>%
        pivot_longer(
                cols = c(Shale,Sand,Sandstone,Gravel,Clay,Limestone,Unknown),
                names_to = "Lithology",
                values_to = "Presence"
        )

# PLOTTING

plot_dp_relabund_nitrate <- ggplot(dp_df2, 
                                aes(y=Nitrate, x=reorder(Site, Nitrate),
                                size = Denitrifier_Rel_Abund)
                                ) + 
            geom_point() +
            geom_hline(yintercept = 0.02,
                linetype="dashed",
                color = "grey40",
                linewidth = 0.6
            )+
            scale_size_continuous(range = c(0.1,14))+
            ylab("Denitrifier Relative\nAbundance vs Nitrate\n\nNitrate (mg/L)") +
            labs(size = "Denitrifier MAGs\nRelative Abundance")+
            guides(
                size = guide_legend(
                        title.position = "top",
                        direction = "horizontal"
                )
            ) +
            scale_y_log10(
                labels = label_number(),
                limits = c(0.005,NA)
            ) +
                theme_minimal() +
                theme(
                        axis.title.x = element_blank(),
                        axis.title.y = element_text(size = 15),
                        axis.ticks = element_blank(),
                        axis.line = element_blank(),
                        panel.grid = element_blank(),
                        axis.text.y = element_text(size = 15),
                        legend.position = c(0.2,0.7),
                        legend.background = element_rect(fill="white",color="black"),
                        legend.text = element_text(size = 10),
                        legend.key.size = unit(0.5,'cm')
                )

plot_dp_relabund_nitrate

#ggsave("plot_dp_relabund_nitrate.png", plot = plot_dp_relabund_nitrate, width = 6, height = 6, dpi = 300)

plot_dp_index_nitrate <- ggplot(dp_df2, 
                                aes(y=Nitrate, x=reorder(Site, Nitrate),
                                size = Index)
                                ) + 
            geom_point() +
            scale_y_log10(
                labels = c("0", "0.01", "0.1", "1", "10", "100"),  # include "0" label
                breaks = c(0.001, 0.01, 0.1, 1, 10, 100),           # fake “0” at ~0.005
                limits = c(0.001,NA)
            )+
            geom_hline(yintercept = 0.02,
                linetype="dashed",
                color = "grey40",
                linewidth = 0.6
            )+
            scale_size_continuous(range = c(0.1,14))+
            ylab("Denitrification Index\nvs Nitrate\n\nNitrate (mg/L)") +
            labs(size = "% Denitrifier MAGs / Total MAGs Abundance")+
            guides(
                size = guide_legend(
                        title.position = "top",
                        direction = "horizontal"
                )
            ) +
            scale_y_log10(
                labels = label_number(),
                limits = c(0.005,NA)
            ) +
                theme_minimal() +
                theme(
                        axis.title.x = element_blank(),
                        axis.title.y = element_text(size = 15),
                        axis.ticks = element_blank(),
                        axis.line = element_blank(),
                        panel.grid = element_blank(),
                        axis.text.y = element_text(size = 15),
                        legend.position = c(0.2,0.7),
                        legend.background = element_rect(fill="white",color="black"),
                        legend.text = element_text(size = 10),
                        legend.key.size = unit(0.5,'cm')
                )

plot_dp_index_nitrate

#ggsave("plot_dp_index_nitrate.png", plot = plot_dp_index_nitrate, width = 6, height = 6, dpi = 300)

plot_dp_lithology <- ggplot(
                        dp_df3, aes(y=Lithology, x=reorder(Site,Nitrate))
                        ) +
                geom_point(aes(size=Presence*100, fill = Lithology),
                        shape = 21, 
                        color = "black"
                        ) +
                scale_fill_manual(values = c(
                        'Coal'     = "#0073C2FF",
                        'Sand'      = "#EFC000FF",
                        'Sandstone' = "#868686FF",
                        'Shale'    = "#CD534CFF",
                        'Gravel'     = "#008B8BFF",
                        'Clay'      = "#7E6148FF",
                        'Limestone' = "#6DBE45FF",
                        'Unknown' = "grey80"
                ) ) +
                ylab("Aquifer Lithology\n\n\n") +
                #scale_y_discrete(expand = c(2, 0) ) +
                theme_minimal() +
                theme(
                        axis.title.x = element_blank(),
                        axis.title.y = element_text(size =15),
                        axis.ticks = element_blank(),
                        axis.line = element_blank(),
                        panel.grid = element_blank(),
                        axis.text.y = element_text(size = 15),
                        legend.position = "none"#,
                        #plot.margin = margin(0,0,0,0)
                )


# Plot well depth

plot_welldepth <- ggplot(dp_df2,
                        aes(y=Well_Depth, x=reorder(Site,Nitrate))
                ) +
                geom_col(
                        fill="#0073C2FF",
                        color="black",
                        width = 0.6
                ) +
                scale_y_reverse() +
                theme_minimal()+
                ylab("Well Depth\n\n\nmeters(m)") +
                theme(
                        axis.title.x = element_blank(),
                        axis.title.y = element_text(size =15),
                        axis.ticks = element_blank(),
                        axis.line = element_blank(),
                        panel.grid = element_blank(),
                        axis.text.y = element_text(size = 15),
                        legend.position = "none"#,
                        #plot.margin = margin(0,0,0,0)
                )


plot_welldepth

# Plots Combined 

plot_combined_1 <- plot_dp_relabund_nitrate / plot_spacer() / plot_welldepth / plot_spacer() /plot_dp_lithology +
        plot_layout(heights = unit(c(5,0.5,3,0.5,4.5), c('cm','cm','cm','cm','cm')),
                #guides = "collect",
                axes = "collect"
        ) +
        plot_annotation(
                tag_levels = "A") &
        theme(
                plot.tag = element_text(size = 20),
                plot.margin = margin(1,0,0,1),
                panel.spacing = unit(0, "lines"),
                axis.text.x = element_text(
                        size = 15,
                        angle = 90,
                        hjust = 1
                )
        )

plot_combined_1

plot_combined_2 <- plot_dp_index_nitrate / plot_spacer() / plot_welldepth / plot_spacer() /plot_dp_lithology +
        plot_layout(heights = unit(c(5,0.5,3,0.5,4.5), c('cm','cm','cm','cm','cm')),
                #guides = "collect",
                axes = "collect"
        ) +
        plot_annotation(
                tag_levels = "A") &
        theme(
                plot.tag = element_text(size = 20),
                plot.margin = margin(1,0,0,1),
                panel.spacing = unit(0, "lines"),
                axis.text.x = element_text(
                        size = 15,
                        angle = 90,
                        hjust = 1
                )
        )

plot_combined_2


#ggsave('combined_plot_denitrifier_relabund.png', plot = plot_combined_1, height = 10, dpi = 300)
#ggsave('combined_plot_denitrifier_index.png', plot = plot_combined_2, height = 10, dpi = 300)

