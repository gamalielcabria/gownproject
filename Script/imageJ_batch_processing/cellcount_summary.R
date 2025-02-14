library("tidyverse")
library("ggplot2")
setwd("/home/glbcabria/Workbench/P3/Results/CellcountPhotos/")
getwd()

# Read the CSV file
csv_D1<-read_csv("/home/glbcabria/Workbench/P3/Results/CellcountPhotos/D-1/Edited/D-1.summary.csv", col_names=TRUE)
csv_D2<-read_csv("/home/glbcabria/Workbench/P3/Results/CellcountPhotos/D2/Edited/D2.summary.csv", col_names=TRUE)
csv_D8<-read_csv("/home/glbcabria/Workbench/P3/Results/CellcountPhotos/D8/Edited/D8.summary.csv", col_names=TRUE)
csv_D14<-read_csv("/home/glbcabria/Workbench/P3/Results/CellcountPhotos/D14/Edited/D14.summary.csv", col_names=TRUE)

csv_D1<-csv_D1 %>% 
    mutate(Slice=gsub("SMPL","D0-SMPL",Slice)) %>%
    mutate(DF=0.5)

csv_D2<-csv_D2 %>% mutate(DF=0.5)
csv_D8<-csv_D8 %>% mutate(DF=0.5) 
csv_D14<-csv_D14 %>% mutate(DF=1.0)

cellcountdf<-rbind(csv_D1,csv_D2,csv_D8,csv_D14)
cellcountdf<-cellcountdf %>%
    filter(Slice != "Slice" & Use != "N") %>%
    mutate(CorrectCount = ifelse( is.na(Corrected), Count, Corrected)) %>%
    mutate(FinalCount = ( as.numeric(CorrectCount) * DF) / (6.25e-6)) %>%
    separate(Slice, c('Day', 'SMPL', 'DAPI', 'Experiment', 'Mag', 'Box'), sep = '-') %>%
    select(-c("Total Area","Average Size","%Area",
        "Mean","Corrected","Count","DF","CorrectCount","Use",
        "SMPL", "DAPI", "Mag"
        )
    ) 

count_summary<-cellcountdf %>%
    mutate(Experiment = case_when(
        Day == "D8" & Experiment == "GR1" ~ "OR1",
        Day == "D8" & Experiment == "TGR1" ~ "GR1",
        TRUE ~ Experiment
    )) %>%
    mutate(Rock.Type = case_when(
        str_detect(Experiment, "G") ~ "Sand",
        str_detect(Experiment, "Y") ~ "Shale",
        str_detect(Experiment, "O") ~ "Sandstone",
        str_detect(Experiment, "B") ~ "Coal",
        TRUE ~ "other"
    )) %>%
    mutate(Setup = case_when(
        str_detect(Experiment, "P") ~ "LN-NO",
        str_detect(Experiment, "W") ~ "NN-NO",
        str_detect(Experiment, "R") ~ "LN-WO",
        TRUE ~ "other"
    )) %>%
    group_by(Day, Rock.Type, Setup) %>%
    summarise(
        Mean.Count=mean(FinalCount),
        SD.Count=sd(FinalCount),
        NumBoxes=n()
    ) %>%
    mutate(
        Day = as.numeric(gsub("D",'',Day)),
    )

# Volume of the 1/16 of the corner cells:
# 0.25mm*0.25*0.10mm = 0.00625 cubic mm = 6.25e-6 ml

plot_count<-ggplot(count_summary) +
    geom_line(aes(x=Day, y=Mean.Count, color=Setup)) +
    geom_point(aes(x=Day, y=Mean.Count, color=Setup), size=2) +
    facet_wrap(~Rock.Type) +
    scale_color_discrete(
        labels = c(
        "LN-NO" = "0.5mM Nitrate - No Oxygen",
        "LN-WO" = "0.5mM Nitrate - With Oxygen",
        "NN-NO" = "No Nitrate - No Oxygen") 
        )
plot_count

ggsave(plot = plot_count, filename = "/home/glbcabria/Workbench/P3/Results/plot_count_D1-D14.png",
    width = 2000,
    height = 2000,
    units = "px",
    dpi = 300
    )
