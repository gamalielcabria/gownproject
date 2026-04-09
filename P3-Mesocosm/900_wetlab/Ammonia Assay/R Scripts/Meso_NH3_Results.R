#Meso Experiment Ammonia Assay results
#March 13th 2026
#Madison Laing

library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(readr)

setwd("/home/gam/github/gownproject/P3-Mesocosm/900_wetlab/Ammonia Assay")

NH3_data <- read_csv("Tables/Ammonia_Assay.csv", col_names = TRUE)



#Calculating mean and standard deviation
NH3_data_summary <- NH3_data %>%
  group_by(Day, Rock_Type) %>%
  summarise(
    Mean_NH3 = mean(Amount_NH3, na.rm = TRUE),
    SD_NH3 = sd(Amount_NH3, na.rm = TRUE),
    n = n()
  ) %>%
  mutate(Day = as.numeric(Day),
         Mean_NH3 = as.numeric(Mean_NH3),
         Chemical = as.character("Ammonium"))%>%
  ungroup()

#Plotting
Meso_NH3_Plot <- ggplot(NH3_data_summary, aes(x = Day, y = Mean_NH3, color = Rock_Type)) +
  geom_point(size = 3) +
  geom_line(aes(group = Rock_Type), linewidth = 1) +
  geom_errorbar(aes(ymin = Mean_NH3 - SD_NH3, ymax = Mean_NH3 + SD_NH3), width = 0.2) +
  facet_wrap(~ Chemical, scales = "free_y") +
  #facet_wrap(~ Rock_Type, scales = "fixed") +
  scale_color_manual(
    values = c(
      "Coal"      = "#0073C2",
      "Sand"      = "#EFC000",
      "Sandstone" = "#868686",
      "Shale"     = "#CD534C",
      "No_Rocks"  = "black"
    )
  ) +
  #scale_x_continuous(breaks = 0:16, limits = c(0, 16))+
  theme_minimal() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 18),
    # legend.position = c(1, 0),
    # legend.justification = c(1, 0),
    # legend.text = element_text(size = 18),
    # legend.title = element_blank()
    legend.position = "none"
  )+  
  labs(
    x = "Day",
    y = "Amount of NH3 (mg/L)",
  )
Meso_NH3_Plot

ggsave(filename = "Meso_NH3_Plot.png", plot = Meso_NH3_Plot, width = 6, height = 6, dpi = 300, bg = "white")
#write_csv(pH_data_summary, "C:/Users/madis/GOWN_LeachingExp/pH Results/pH_data_summary.csv")
