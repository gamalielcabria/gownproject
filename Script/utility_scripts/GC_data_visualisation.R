library("tidyverse")
library(scales)
library(ggpubr)
library(rstatix)
library(patchwork)
library(stringr)
library(ggsci)

# args <- commandArgs(trailingOnly = TRUE)
# input<-as.character(args[1])
input <- '/home/glbcabria/Workbench/P3/Results/GCPhotos/D2_GC.csv'
input_ssa <- '/home/glbcabria/Workbench/P3/Results/GCPhotos/SurfaceAreaPerReplicate.csv'
input_ccs <- '/home/glbcabria/Workbench/P3/Results/count_summary.csv'
dataframe <- read_csv(input, col_names = TRUE)

gc_df<- dataframe %>%
  filter(ImageName != "ImageName") %>%
  separate(Filename, into = c("N2O","Date_GC", "Time_GC","SampleName"), sep="_") %>%
  #filter(N2O != "CALIBRATION") %>%
  filter(grepl("SAMPL|SMPL|SAMPLE|SM|SMPI|CTRL|TRL",SampleName)) %>%
  separate(SampleName, into = c("Sample", "Date_Start", "Time_Start", "Type"), sep="-") %>%
  mutate(Time_GC=gsub("(^[0-9]{3})(A|P)M","0\\1 \\2M",Time_GC, fixed = FALSE, perl =TRUE),
         Time_Start=gsub("(^[0-9]{3})(A|P)M","0\\1 \\2M",Time_Start, fixed = FALSE, perl =TRUE),
         Time_GC=gsub("(^[0-9]{4})(A|P)M","\\1 \\2M",Time_GC, fixed = FALSE, perl =TRUE),
         Time_Start=gsub("(^[0-9]{4})(A|P)PM","\\1 \\2M",Time_Start, fixed = FALSE, perl =TRUE)
         ) %>%
  unite("DateTime_GC", Date_GC:Time_GC, sep = ' ') %>%
  mutate(DateTime_GC = strptime(DateTime_GC, "%Y%m%d %I%M %p")) %>%
  unite("DateTime_Start", Date_Start:Time_Start, sep = ' ') %>%
  mutate(DateTime_Start = strptime(DateTime_Start, "%Y%m%d %I%M %p")) %>%
  mutate(IncTime = DateTime_GC-DateTime_Start) %>%
  mutate(Amount = ifelse(grepl("\\d\\.\\d", Amount, perl=TRUE), Amount, 0)) %>%
  mutate(Amount = as.numeric(Amount)) %>%
  mutate(Type =  gsub("^SAMPL$|^SM$|^SMPI$", "SMPL", Type),
         Type = gsub("^TRL$|^TR$|^TRI$|^CTR$", "CTRL", Type)) %>%
  select(Amount, Sample, Type, IncTime) %>% #
  mutate(Rock = case_when(
    grepl("^CO", Sample) ~ "Coal",
    grepl("^SH", Sample) ~ "Shale",
    grepl("^SS|^SN", Sample) ~ "Sandstone",
    grepl("^SA", Sample) ~ "Sand",
    .default = "CTRL"
  )) %>%
  mutate(Experiment = case_when(
    grepl("LNWO\\d$|LNW0\\d$", Sample) ~ "Low Nitrate-With Oxygen",
    grepl("LNNO\\d$|LNN0\\d$", Sample) ~ "Low Nitrate-No Oxygen",
    grepl("NNWO\\d$|NNW0\\d$", Sample) ~ "No Nitrate-With Oxygen",
    grepl("NNNO\\d$|^[A-Z]{2}NNO\\d$|NNN0\\d$|^[A-Z]{2}NN0\\d$", Sample) ~ "No Nitrate-No Oxygen",
    .default = "CTRL"
  )) %>%
  mutate(Replicate = str_extract(Sample, "(\\d)$")) %>%
  unite("SampleID",Rock:Replicate, sep = "-",remove = FALSE) %>%
  select(-Sample)

####################
# Rate Calculation #
####################
# Compute the individual slopes
gc_df_indiv_slopes <- gc_df %>%
  group_by(SampleID, Rock, Experiment) %>%
  summarise(slope = coef(lm(Amount ~ IncTime))[2],  # Extract slope
            intercept = coef(lm(Amount ~ IncTime))[1],  # Extract intercept
            .groups = "drop")

# Compute the aggregate slopes
gc_df_aggregate_slopes <- gc_df %>%#filter(Amount != 0)%>%
  filter(!grepl('CTRL',SampleID))%>%
  group_by(Rock, Experiment) %>%
  summarise(slope = coef(lm(Amount ~ IncTime))[2],  # Extract slope
            intercept = coef(lm(Amount ~ IncTime))[1],  # Extract intercept
            .groups = "drop")

# Normalize by Surface Area
ssa_df <- read_csv(input_ssa, col_names = TRUE)
ssa_df2 <- ssa_df %>% 
  mutate(`Sample Name` = gsub('-','',`Sample Name`) ) %>%
  rename(Sample = `Sample Name`, Surface.Area = `Surface Area cm2`) %>%
  mutate(Rock = case_when(
    grepl("^CO", Sample) ~ "Coal",
    grepl("^SH", Sample) ~ "Shale",
    grepl("^SS|^SN", Sample) ~ "Sandstone",
    grepl("^SA", Sample) ~ "Sand",
    .default = "CTRL"
  )) %>%
  mutate(Experiment = case_when(
    grepl("LNWO\\d$|LNW0\\d$", Sample) ~ "Low Nitrate-With Oxygen",
    grepl("LNNO\\d$|LNN0\\d$", Sample) ~ "Low Nitrate-No Oxygen",
    grepl("NNWO\\d$|NNW0\\d$", Sample) ~ "No Nitrate-With Oxygen",
    grepl("NNNO\\d$|^[A-Z]{2}NNO\\d$|NNN0\\d$|^[A-Z]{2}NN0\\d$", Sample) ~ "No Nitrate-No Oxygen",
    .default = "CTRL"
  )) %>%  
  mutate(Replicate = str_extract(Sample, "(\\d)$")) %>%
  unite("SampleID", Rock:Replicate, sep = '-',remove = TRUE,)
  
normalized_gc_df_slopes <- gc_df_indiv_slopes %>%
  left_join(ssa_df2 %>% select(SampleID, Surface.Area), by = 'SampleID') %>%
  filter(!grepl('CTRL',SampleID)) %>%
  mutate(Normalized.Rate = slope / Surface.Area) 

normalized_gc_df_agg_slopes <- normalized_gc_df_slopes %>%
  group_by(Rock, Experiment) %>%
  summarise(Mean.Rate.N2OperSA = mean(Normalized.Rate),
            SD.Rate.N2OperSA = sd(Normalized.Rate),
            Mean.Intercept = mean(intercept))

# Normalized by Cell count (Replicate 1 only)
ccs_df <- read_csv(input_ccs, col_names = TRUE)
ccs_df2 <- ccs_df %>% 
  mutate(Setup = gsub('-','',Setup)) %>%
  mutate(Experiment = case_when(
    grepl("LNWO", Setup) ~ "Low Nitrate-With Oxygen",
    grepl("LNNO", Setup) ~ "Low Nitrate-No Oxygen",
    grepl("NNWO", Setup) ~ "No Nitrate-With Oxygen",
    grepl("NNNO", Setup) ~ "No Nitrate-No Oxygen",
    .default = "CTRL"
  )) %>%  
  mutate(Replicate = 1) %>%
  unite("SampleID", c('Rock.Type','Experiment','Replicate'), sep = '-',remove = TRUE,) %>%
  filter(Day == 2) %>% 
  inner_join(gc_df_indiv_slopes, by = 'SampleID')

normalized_gc_df_slopes_ccs <- gc_df_indiv_slopes %>%
  left_join(ccs_df2 %>% select(SampleID, Mean.Count), by = 'SampleID') %>%
  filter(!grepl('CTRL',SampleID)) %>%
  mutate(Normalized.Rate = slope / Mean.Count) %>%
  filter(!is.na(Normalized.Rate))

#################
# Visualization #
#################

# Setup
# Total Moles in Gas  PV=nRT
# R = universal gas constant (0.0821 L·atm/(mol·K))
temp = 22 #celsius
press = 1.5 #atm
volume = 0.110 # L
moles_total = (press * volume) / (0.0821 * (22+273.15))

# Plot Normalized Slopes
plot_DNP_agg_norm  <- ggplot(gc_df%>%#filter(Amount != 0)%>%
                          filter(!grepl('CTRL',SampleID)), 
                        aes(x=IncTime, y=Amount))+
  geom_point(aes(color=Rock, shape=Experiment), size=3)+
  geom_abline(data = gc_df_aggregate_slopes, aes(slope = slope, intercept = intercept), 
              color = "black", linetype = "dashed", linewidth = 1) +
  facet_grid(Rock ~ Experiment)+
  scale_y_continuous(
    labels = function(x) label_comma() ((x/100) * moles_total * 1e9 ), 
    name = "nmol N2O"
  ) +
  scale_color_jco() +
  xlab("Incubation Time (hrs)") +
  theme_light() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.position = "none"#,
    #aspect.ratio = 850/790
    )
plot_DNP_agg_norm


###################
# Plot Statistics #
###################

# Conversion of Mole % to mmol: mmol = (mole%/100) * moles_total * 1000

# Normalized rate
  stat_gc_norm <- normalized_gc_df_slopes %>%
    mutate(Normalized.Rate = (Normalized.Rate/100) * moles_total * 1e9 ) %>% 
    mutate(xint=(0-intercept)/slope) %>%
    filter(!grepl('CTRL',SampleID))

# Stat of Normalized Slopes
  stat_plot_norm <- ggboxplot(stat_gc_norm,
                         x = 'Rock',
                         y = 'Normalized.Rate',
                         color = 'Rock',
                         add = c("jitter", "mean_sd"),
                         palette = "jco") +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +    
    stat_compare_means(comparisons = list(c("Coal", "Sand"), c("Shale", "Coal"), c("Sandstone", "Coal"),
                                          c("Sandstone", "Sand"), c("Shale", "Sand"),c("Shale", "Sandstone")
                                          ), 
                       method = "wilcox.test", 
                       paired = FALSE, #hide.ns = TRUE,
                       label="p.signif")+
    ylab("Denitrification Rate\n(nmol N2O per hour per cm^2)" ) +
    theme_minimal()+
    theme(
      axis.text = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.title.x = element_text(size = 18) 
    )
  stat_plot_norm

# Normalized by cell count (replicate 1)
  stat_gc_norm_cc <- normalized_gc_df_slopes_ccs %>%
    mutate(Normalized.Rate = (Normalized.Rate/100) * moles_total * 1e9 ) %>% 
    mutate(xint=(0-intercept)/slope) %>%
    filter(!grepl('CTRL',SampleID))

  # Stat of Normalized Slopes by Cell Count
  stat_plot_norm_cc <- ggboxplot(stat_gc_norm_cc,
                              x = 'Rock',
                              y = 'Normalized.Rate',
                              color = 'Rock',
                              add = c("jitter", "mean_sd"),
                              palette = "jco") +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
    stat_compare_means(comparisons = list(c("Coal", "Sand"), c("Shale", "Coal"), c("Sandstone", "Coal"),
                                          c("Sandstone", "Sand"), c("Shale", "Sand"),c("Shale", "Sandstone")
                                          ), 
                       method = "wilcox.test", 
                       paired = FALSE, #hide.ns = TRUE,
                       label="p.signif") +
    ylab("Denitrification Rate\n(nmol N2O per hour per cell)" ) +
    theme_minimal()+
    theme(
      axis.text = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.title.x = element_text(size = 18) 
    )
  stat_plot_norm_cc
  
# Aggregate
  layout <- "
  AAABB
  AAABB
  "
  plot_norm_cc_sa_rate <- plot_DNP_agg_norm + ( stat_plot_norm / stat_plot_norm_cc) + 
    plot_annotation(tag_levels = "A") +
    plot_layout(guides = "collect", design = layout) & 
    theme(
      legend.position = "none",
      plot.tag = element_text(size = 20)
      )

  # plot_norm_cc_sa_rate <- ggarrange( plot_DNP_agg_norm, stat_plot_norm, stat_plot_norm_cc, ncol=3, nrow=1,
  #                                    common.legend = TRUE, legend = "right",
  #                                    labels = c("A","B","C"),
  #                                    font.label = list(size =14, face = "bold"))
  plot_norm_cc_sa_rate
  
  ggsave("/home/glbcabria/Workbench/P3/Results/GCPhotos/plot_norm_cc_sa_rate.png", plot_norm_cc_sa_rate,
         width = 3900, height = 2600,
         units = c('px'),
         dpi = 200
         )

  
#=================================================================#
# # Plot individual slopes
# plot_DNP_indiv <- ggplot(gc_df%>%filter(Amount != 0), 
#                    aes(x=IncTime, y=Amount))+
#   geom_point(aes(color=Rock, shape=Experiment))+
#   geom_abline(data = gc_df_indiv_slopes, aes(slope = slope, intercept = intercept), 
#               color = "black", linetype = "dashed", linewidth = 1) +
#   #facet_grid(Rock ~ Experiment, scales = "free")
#   facet_wrap(~SampleID, scales = "free")
# plot_DNP_indiv
# 
# # Plot Aggregate slopes
# plot_DNP_agg  <- ggplot(gc_df%>%#filter(Amount != 0)%>%
#                           filter(!grepl('CTRL',SampleID)), 
#                          aes(x=IncTime, y=Amount))+
#   geom_point(aes(color=Rock, shape=Experiment))+
#   geom_abline(data = gc_df_aggregate_slopes, aes(slope = slope, intercept = intercept), 
#               color = "black", linetype = "dashed", linewidth = 1) +
#   facet_grid(Rock ~ Experiment)+
#   theme_light() +
#   ylab("%Mole N2O")+
#   theme(
#     axis.text = element_text(size = 15)
#   )
# plot_DNP_agg
# # Stat of non-normalized
# stat_gc <- gc_df_indiv_slopes %>% 
#   mutate(slope = slope * 1000) %>% 
#   mutate(xint=(0-intercept)/slope) %>%
#   filter(!grepl('CTRL',SampleID))
# 
# stat_plot <- ggboxplot(stat_gc,
#                        x = 'Rock',
#                        y = 'slope',
#                        color = 'Rock',
#                        add = c("jitter", "mean_sd"),
#                        palette = "jco") +
#   stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
#   stat_compare_means(comparisons = list(c("Coal", "Sand"), c("Sandstone", "Sand"), 
#                                         c("Shale", "Sand"),c("Shale", "Sandstone")), 
#                      method = "wilcox.test",
#                      label="p.signif") +
#   # stat_compare_means(method = "anova", label.y = max(df$Value) + 1) +
#   ylab("Denitrification Rate (ppm N2O per hour)") +
#   theme_minimal()
#   
# stat_plot
# 
# stat_plot2 <- ggboxplot(stat_gc,
#                        x = 'Experiment',
#                        y = 'slope',
#                        color = 'Experiment',
#                        add = c("jitter", "mean_sd"),
#                        palette = "jco") +
#   stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
#   stat_compare_means(comparisons = list(c("Low Nitrate-With Oxygen", "Low Nitrate-No Oxygen"), c("No Nitrate-No Oxygen", "Low Nitrate-No Oxygen") ), 
#                      method = "wilcox.test",
#                      label="p.signif") +
#   # stat_compare_means(method = "anova", label.y = max(df$Value) + 1) +
#   ylab("Denitrification Rate (ppm N2O per hour)") +
#   theme_minimal() + 
#   facet_grid(~Rock)
# stat_plot2

# # Output the table
# gc_processed_output <- paste0("Table_Processed", basename(input))
# write_csv(gc_df, file = gc_processed_output)
# 
# gc_slopes_indiv_output <- paste0("Table_Slopes", basename(input))
# write_csv(gc_df_indiv_slopes, file = gc_slopes_indiv_output)
# 
# gc_slopes_aggregate_output <- paste0("Table_Slopes_Aggregate", basename(input))
# write_csv(gc_df_aggregate_slopes, file = gc_slopes_aggregate_output)
# 
# output_plot1 <-paste0("Plot_DNP_indiv_", str_split_1(basename(input), '\\.')[1],'.png')
# ggsave(output_plot1, plot_DNP_indiv, width = 4000, height = 3000, units = "px",dpi=300)
# 
# output_plot2 <-paste0("Plot_DNP_aggregte_", str_split_1(basename(input), '\\.')[1],'.png')
# ggsave(output_plot2, plot_DNP_agg, width = 2000, height = 3000, units = "px",dpi=300)