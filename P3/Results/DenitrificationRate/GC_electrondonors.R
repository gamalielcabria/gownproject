## Libraries
#install.packages("pdftools")
library(pdftools)
library(tidyverse)
library(broom)
library(gridExtra)
library(patchwork)
library(minpack.lm)

##################### Functions #############################################

find_N2O_firstmatch <- function(text, lower = 1.19, upper = 1.3) {
  # Split on literal "\n"
  txt_lines <- strsplit(text, "\\\n")[[1]]
  
  # Extract numbers from each line
  all_numbers <- regmatches(txt_lines, gregexpr("\\d+\\.\\d+", txt_lines))
  
  for (i in seq_along(all_numbers)) {
    # Only check lines containing "N2O" (case-insensitive)
    if (grepl("N2O|?", txt_lines[i], ignore.case = TRUE)) {
      nums <- as.numeric(all_numbers[[i]])
      match_nums <- nums[nums >= lower & nums <= upper]
      
      if (length(match_nums) > 0) {
        return(txt_lines[i])  # return the full line
      }
    }
  }
  
  return(NULL)
}

############################# PROCESSING INPUTS #############################

# Working Directory
input <- '/home/gam/github/gownproject/P3/Results/DenitrificationRate/combined'
#input_day <- 'D87LIQ' # 'FECH4MNH2'

# Determine if Amount is in PPM or PCT MOLE
amount_unit <- 'PPM' #'PCT'

# Set the Retention Time for N2O
lower_bound <- 1.79
upper_bound <- 1.9

setwd(input)

# File Location
files <- list.files(pattern="\\.pdf$", ignore.case = TRUE)

# Create empty dataframe to store results
GC_df <- data.frame(File = character(),
                         Results = character(),
                         stringsAsFactors = FALSE)

# Read through files
for (file in files) {
  cat(file)
  txt <- pdf_text(file)
  #filename <- regmatches(txt[2], regexpr("\\d{8}[-_][A-Za-z0-9_-]+\\.D", txt[2]))
  filename <- gsub("\\.pdf",'',as.character(file))
  N2O_line <- trimws(find_N2O_firstmatch(txt[2], lower = lower_bound, upper = upper_bound))
  
  results <- gsub("\\s+", "_", N2O_line)
  #cat(results)
  GC_df <- rbind( GC_df,
                  data.frame(File = filename,
                             Results = results,
                             stringsAsFactors =  FALSE)
    
  )
}

GC_processed_df <- GC_df %>%
  mutate(Amount = ifelse(grepl('POSTPROCESS',File,ignore.case = TRUE), NA,
                         as.numeric(map_chr(str_split(Results, "_"), ~ .x[length(.x) -1]) )
                         )
         ) %>%
  mutate(Area = ifelse(grepl('POSTPROCESS',File,ignore.case = TRUE),
                       as.numeric(map_chr(str_split(Results, "_"), ~ .x[length(.x) -2]) ),
                       as.numeric(map_chr(str_split(Results, "_"), ~ .x[length(.x) -3]) ) 
                       ) 
         )%>%
  filter(!grepl("PREPROCESS",File, ignore.case = TRUE)) %>%
  mutate(File = gsub("-POSTPROCESS", "", File)) %>%
  separate(File, into = c('Date-Time','Type','Sample'), sep = "_") %>%
  select(-Results) %>%
  separate('Date-Time', into = c("Date","Time"), sep = '-' ) %>%
  #mutate(Amount = ifelse(Type != 'Calibration', NA, Amount) ) 
  mutate(Amount = ifelse(grepl('Standards-0', Sample, ignore.case = TRUE), 0, Amount ),
         Area = ifelse(grepl('Standards-0', Sample, ignore.case = TRUE), 0, Area ) 
         ) %>%
  mutate(Sample = gsub("-\\d+$|-\\d{1}[A-Za-z]{2}$|-FILTERED", "", Sample) ) %>% # Remove the -8PM/-1/-2 from Sample Name
  replace_na(list(Area = 0))

# Checking the Amount whether in PPM or PCT MOLE
# Convert it to mmol/L N2O
if (amount_unit == "PCT") {
  GC_processed_df <- GC_processed_df %>%
    mutate(Amount = Amount * 10000 * 0.0227206)
} else if (amount_unit == "PPM") {
    GC_processed_df <- GC_processed_df %>%
        mutate(Amount = Amount * 0.02272066) # Convert PPM to mmol/L N2O
}
###################### Process Calibration Models ###########################

# Creating Calibration Model
models <- GC_processed_df %>%
  filter(grepl("^Calibration$", Type, ignore.case = TRUE)) %>%
  group_by(Date) %>%
  summarise(
    model = list(lm(Amount ~ Area, data = cur_data())), 
    .groups = "drop"
  )

# Join models back to main dataframe by Date
GC_filled <- GC_processed_df %>%
  filter(grepl("^Sample$", Type, ignore.case = TRUE)) %>%
  left_join(models, by = "Date") %>%
  mutate(
    Amount = map2_dbl(model, Area, ~ {
      if (is.null(.x) || all(is.na(coef(.x)))) { 
        NA_real_ # no matching models for that Date
      } else {
        as.numeric(predict(.x, newdata = data.frame(Area = .y)))
      } 
    })
  ) %>%
  select(-model) %>%
  separate(Sample, into = c('Start_Date', 'Start_Time', 'EXPT'), sep = '-',extra = 'merge', fill = 'right') %>%
  mutate(
    # Setup: detect starting letters in EXPT
    Setup = case_when(
      str_detect(EXPT, "SS")    ~ "Sandstone",
      str_detect(EXPT, "SH")    ~ "Shale",
      str_detect(EXPT, "SA")    ~ "Sand",
      str_detect(EXPT, "CO")    ~ "Coal",
      str_detect(EXPT, "CL")    ~ "Control",
      str_detect(EXPT, "BLANK|MQ")  ~ "Blank",
      TRUE                       ~ NA_character_
    ),
    Setup2 = case_when(
      str_detect(EXPT, "ORG")    ~ "Sandstone",
      str_detect(EXPT, "INO")    ~ "Shale",
      str_detect(EXPT, "CTL")    ~ "Control",
      TRUE                       ~ "Blank"
    ),
    # Replicate: last character of EXPT
    Replicate = str_extract(EXPT, "[0-9B]$"),
    Replicate = if_else(Setup == "Blank" & is.na(Replicate), "B", Replicate)
   ) 

######################## Processing ROCK COMPARISON #########################
# Filter SED only data
GC_sed_df <- GC_filled %>%
  filter(!is.na(EXPT) & str_detect(EXPT, regex('D88SED', ignore_case = TRUE))) %>%
  mutate(
    EndDT   = parse_date_time(paste(Date, '-', Time), orders = "%Y%m%d - %I%M%p", tz = "UTC"),
    StartDT = parse_date_time(paste(Start_Date, '-', Start_Time), orders = "%Y%m%d - %I%M%p", tz = "UTC"),
    Total_Time = as.numeric(difftime(EndDT, StartDT, units = "hours"))
  )


########## Final Table for Publication ##############
Table_GC_Sediments <- GC_sed_df %>%
  select(Total_Time,Amount,Setup,Replicate) %>%
  as.data.frame()
  

######################### Adding Best Fit Line ##############################

# 1. Filter data (remove 'B' replicates and 'Blank' Setup)
filtered_data <- Table_GC_Sediments %>%
  filter(!is.na(Replicate)) %>%
  filter(Replicate != "B", Setup != "Blank") %>%
  filter(Replicate != 4)  
  
# 2. Define fitting function (same as before)
fit_models <- function(df) {
  start_linear <- list(a = min(df$Amount), b = 0)
  fit_linear <- try(lm(Amount ~ Total_Time, data = df), silent = TRUE)  
  aic_linear <- if(inherits(fit_linear, "lm")) AIC(fit_linear) else NA
  aic_vals <- c(linear = aic_linear)
  best_model <- names(which.min(aic_vals))
  list(
    fit_linear = fit_linear,
    best_model = best_model,
    aic = aic_vals
  )
}

# 3. Fit models by Setup and keep the best fit object
fit_results <- filtered_data %>%
  group_by(Setup) %>%
  nest() %>%
  mutate(fits = map(data, fit_models)) 

# 4. Create prediction data frames with best fit per Setup
predictions <- fit_results %>%
  mutate(
    pred_df = pmap(list(data, fits, Setup), function(data, fits, Setup) {
      # Define a sequence for Total_Time over observed range for smooth curve
      new_time <- seq(min(data$Total_Time), max(data$Total_Time), length.out = 100)
      
      # Predict based on best model
      pred_amount <- switch(fits$best_model,
                            linear = predict(fits$fit_linear, newdata = data.frame(Total_Time = new_time))
                            )
      
      tibble(
        Total_Time = new_time,
        Amount = pred_amount
      )
    })
  ) %>%
  select(Setup, pred_df) %>%
  unnest(pred_df)

#################################### PLot for Pub ############################################

#1. Set Colors to certain levels
setup_cols <- c(
  Coal      = "#0073C2",
  Sand      = "#EFC000",
  Sandstone = "#868686",
  Shale     = "#CD534C"
)

## Matches the levels of colors to the levels of Setup
lvls <- names(setup_cols)
filtered_data$Setup  <- factor(filtered_data$Setup,  levels = lvls)
predictions$Setup    <- factor(predictions$Setup,    levels = lvls)

#2. Generate Main Plot
Plot_DNP <- ggplot() +
  geom_point(data = filtered_data, aes(x = Total_Time, y = Amount, color = Setup)) +
  geom_line(data = predictions, aes(x = Total_Time, y = Amount, color = Setup), linewidth = 1) +
  facet_wrap(~Setup) +
  scale_color_manual(values = setup_cols, limits = lvls) +  # enforce mapping & order
    labs(#title = paste0("EXPT ",input_day),
       x = "Total Time",
       y = "Amount (mmol/L N2O)") +
  theme_minimal() +
  theme(
    aspect.ratio = 1,
    #plot.title = element_text(size = 20, hjust = 0.5),
    axis.text = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.placement = "outside",
    strip.background = element_rect(fill="grey90"),
    strip.text = element_text(face = "bold", size = 18),
    #strip.switch.pad.wrap = unit(1, "cm"),
    legend.position = "none"
  )

Plot_DNP

#3. Compute the Slopes and save into a table
slopes_table <- filtered_data %>%
  group_by(Setup, Replicate) %>%
  summarise(
    # Fit a linear model Amount ~ Total_Time for each group
    slope = if (n() > 1) {
      coef(lm(Amount ~ Total_Time))[2]  # slope coefficient
    } else {
      NA_real_
    },
    r2 = if (n() > 1) {
      summary(lm(Amount ~ Total_Time))$r.squared  # R² value
    } else {
      NA_real_
    },
    .groups = "drop"
  )

################## Table for Pub ################
Table_Slopes <- slopes_table %>%
  ungroup()%>%
  group_by(Setup) %>%
  summarise(
    Mean = round(mean(slope), 3),
    Std_Dev = round(sd(slope), 3),
    r2 = round(mean(r2, na.rm = TRUE), 3),
    Replicates = n()
            ) #%>%
  #mutate(DAY = 'D88SED') # MANUALLY SET THE DAY in line
#################################################


# 4. Alternative output format for plot: ADD INSET TABLE
#tableGrob Themes
minimal_theme <- ttheme_minimal(
  padding = unit(c(6, 10), "pt"),
  core = list(
    fg_params = list(fontsize = 15),
    bg_params = list(fill = NA, col = NA)
  ),
  colhead = list(
    fg_params = list(fontface = "bold", fontsize = 15),
    bg_params = list(fill = NA, col = NA)
  ),
  rowhead = list(
    fg_params = list(fontsize = 15),
    bg_params = list(fill = NA, col = NA)
  ),
)

# 5. Create the table grob with the minimal theme
inset_table <- Table_Slopes %>%
  unite('Slope', Mean, Std_Dev, sep = '±') %>%
  select(Setup, Slope, r2)

table_grob <- tableGrob(
  inset_table, 
  rows = NULL, theme = minimal_theme
  )

Plot_DNP_Table_All <- Plot_DNP + inset_element(
  table_grob,
  #on_top = FALSE,
  align_to = 'panel',
  left = 0.15, bottom = 0.7, right = 0.4, top = 0.7
  #left = 0.75, bottom = 0.10, right = .95, top = .35
)
Plot_DNP_Table_All

########################### END OF ROCK COMPARISON #########################


######################## Processing ELECTRON DONOR #########################
# Filter SED only data
GC_ED_df <- GC_filled %>%
  filter(!is.na(EXPT) & str_detect(EXPT, regex('ELECDON', ignore_case = TRUE))) %>%
  mutate(
    EndDT   = parse_date_time(paste(Date, '-', Time), orders = "%Y%m%d - %I%M%p", tz = "UTC"),
    StartDT = parse_date_time(paste(Start_Date, '-', Start_Time), orders = "%Y%m%d - %I%M%p", tz = "UTC"),
    Total_Time = as.numeric(difftime(EndDT, StartDT, units = "hours"))
  )


########## Final Table for Publication ##############
Table_GC_Sediments <- GC_sed_df %>%
  select(Total_Time,Amount,Setup,Replicate) %>%
  as.data.frame()
  

######################### Adding Best Fit Line ##############################

# 1. Filter data (remove 'B' replicates and 'Blank' Setup)
filtered_data <- Table_GC_Sediments %>%
  filter(!is.na(Replicate)) %>%
  filter(Replicate != "B", Setup != "Blank") %>%
  filter(Replicate != 4)  
  
# 2. Define fitting function (same as before)
fit_models <- function(df) {
  start_linear <- list(a = min(df$Amount), b = 0)
  fit_linear <- try(lm(Amount ~ Total_Time, data = df), silent = TRUE)  
  aic_linear <- if(inherits(fit_linear, "lm")) AIC(fit_linear) else NA
  aic_vals <- c(linear = aic_linear)
  best_model <- names(which.min(aic_vals))
  list(
    fit_linear = fit_linear,
    best_model = best_model,
    aic = aic_vals
  )
}

# 3. Fit models by Setup and keep the best fit object
fit_results <- filtered_data %>%
  group_by(Setup) %>%
  nest() %>%
  mutate(fits = map(data, fit_models)) 

# 4. Create prediction data frames with best fit per Setup
predictions <- fit_results %>%
  mutate(
    pred_df = pmap(list(data, fits, Setup), function(data, fits, Setup) {
      # Define a sequence for Total_Time over observed range for smooth curve
      new_time <- seq(min(data$Total_Time), max(data$Total_Time), length.out = 100)
      
      # Predict based on best model
      pred_amount <- switch(fits$best_model,
                            linear = predict(fits$fit_linear, newdata = data.frame(Total_Time = new_time))
                            )
      
      tibble(
        Total_Time = new_time,
        Amount = pred_amount
      )
    })
  ) %>%
  select(Setup, pred_df) %>%
  unnest(pred_df)

#################################### PLot for Pub ############################################

#1. Set Colors to certain levels
setup_cols <- c(
  Coal      = "#0073C2",
  Sand      = "#EFC000",
  Sandstone = "#868686",
  Shale     = "#CD534C"
)

## Matches the levels of colors to the levels of Setup
lvls <- names(setup_cols)
filtered_data$Setup  <- factor(filtered_data$Setup,  levels = lvls)
predictions$Setup    <- factor(predictions$Setup,    levels = lvls)

#2. Generate Main Plot
Plot_DNP <- ggplot() +
  geom_point(data = filtered_data, aes(x = Total_Time, y = Amount, color = Setup)) +
  geom_line(data = predictions, aes(x = Total_Time, y = Amount, color = Setup), linewidth = 1) +
  facet_wrap(~Setup) +
  scale_color_manual(values = setup_cols, limits = lvls) +  # enforce mapping & order
    labs(#title = paste0("EXPT ",input_day),
       x = "Total Time",
       y = "Amount (mmol/L N2O)") +
  theme_minimal() +
  theme(
    aspect.ratio = 1,
    #plot.title = element_text(size = 20, hjust = 0.5),
    axis.text = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.placement = "outside",
    strip.background = element_rect(fill="grey90"),
    strip.text = element_text(face = "bold", size = 18),
    #strip.switch.pad.wrap = unit(1, "cm"),
    legend.position = "none"
  )

Plot_DNP

#3. Compute the Slopes and save into a table
slopes_table <- filtered_data %>%
  group_by(Setup, Replicate) %>%
  summarise(
    # Fit a linear model Amount ~ Total_Time for each group
    slope = if (n() > 1) {
      coef(lm(Amount ~ Total_Time))[2]  # slope coefficient
    } else {
      NA_real_
    },
    r2 = if (n() > 1) {
      summary(lm(Amount ~ Total_Time))$r.squared  # R² value
    } else {
      NA_real_
    },
    .groups = "drop"
  )

################## Table for Pub ################
Table_Slopes <- slopes_table %>%
  ungroup()%>%
  group_by(Setup) %>%
  summarise(
    Mean = round(mean(slope), 3),
    Std_Dev = round(sd(slope), 3),
    r2 = round(mean(r2, na.rm = TRUE), 3),
    Replicates = n()
            ) #%>%
  #mutate(DAY = 'D88SED') # MANUALLY SET THE DAY in line
#################################################


# 4. Alternative output format for plot: ADD INSET TABLE
#tableGrob Themes
minimal_theme <- ttheme_minimal(
  padding = unit(c(6, 10), "pt"),
  core = list(
    fg_params = list(fontsize = 15),
    bg_params = list(fill = NA, col = NA)
  ),
  colhead = list(
    fg_params = list(fontface = "bold", fontsize = 15),
    bg_params = list(fill = NA, col = NA)
  ),
  rowhead = list(
    fg_params = list(fontsize = 15),
    bg_params = list(fill = NA, col = NA)
  ),
)

# 5. Create the table grob with the minimal theme
inset_table <- Table_Slopes %>%
  unite('Slope', Mean, Std_Dev, sep = '±') %>%
  select(Setup, Slope, r2)

table_grob <- tableGrob(
  inset_table, 
  rows = NULL, theme = minimal_theme
  )

Plot_DNP_Table_ED <- Plot_DNP + inset_element(
  table_grob,
  #on_top = FALSE,
  align_to = 'panel',
  left = 0.15, bottom = 0.7, right = 0.4, top = 0.7
  #left = 0.75, bottom = 0.10, right = .95, top = .35
)
Plot_DNP_Table_ED

###########################