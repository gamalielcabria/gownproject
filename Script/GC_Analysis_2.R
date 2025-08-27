## Libraries
#install.packages("pdftools")
library(pdftools)
library(tidyverse)
library(broom)
library(gridExtra)
library(patchwork)

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

#############################################################################

# Working Directory
input <- '/home/glbcabria/Workbench/P3/Expt2/GC/Train'
input_day <- 'D56'

# Determine if Amount is in PPM or PCT MOLE
amount_unit <- 'PCT'

# Set the Retention Time for N2O
lower_bound <- 1.19
upper_bound <- 1.3

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
if (amount_unit == "PCT") {
  GC_processed_df <- GC_processed_df %>%
    mutate(Amount = Amount * 10000)
}

# Calibration Model
models <- GC_processed_df %>%
  filter(grepl("^Calibration$", Type, ignore.case = TRUE)) %>%
  group_by(Date) %>%
  summarise(
    model = list(lm(Amount ~ Area, data = cur_data()))
  )

# Join models back to main dataframe by Date
GC_filled <- GC_processed_df %>%
  filter(grepl("^Sample$", Type, ignore.case = TRUE)) %>%
  left_join(models, by = "Date") %>%
  mutate(
    Amount = map2_dbl(model, Area, ~ {
      if (is.null(.x)) return(NA_real_) # no matching models for that Date
      predict(.x, newdata = data.frame(Area = .y))
    })
  ) %>%
  select(-model) %>%
  separate(Sample, into = c('Start_Date', 'Start_Time', 'EXPT'), sep = '-',extra = 'merge') %>%
  mutate(
    # Setup: detect starting letters in EXPT
    Setup = case_when(
      str_detect(EXPT, "SS")    ~ "Sandstone",
      str_detect(EXPT, "SH")    ~ "Shale",
      str_detect(EXPT, "SA")    ~ "Sand",
      str_detect(EXPT, "CO")    ~ "Coal",
      str_detect(EXPT, "CL")    ~ "Control",
      str_detect(EXPT, "BLANK")  ~ "Blank",
      TRUE                       ~ NA_character_
    ),
    # Replicate: last character of EXPT
    Replicate = str_extract(EXPT, "[0-9B]$")
  ) %>%
  mutate(
    Replicate = if_else(Setup == "Blank" & is.na(Replicate), "B", Replicate)
  ) %>%
  filter(!is.na(EXPT) & str_detect(EXPT, regex(input_day, ignore_case = TRUE))) %>%
  mutate(
    EndDT   = parse_date_time(paste(Date, '-', Time), orders = "%Y%m%d - %I%M%p", tz = "UTC"),
    StartDT = parse_date_time(paste(Start_Date, '-', Start_Time), orders = "%Y%m%d - %I%M%p", tz = "UTC"),
    Total_Time = as.numeric(difftime(EndDT, StartDT, units = "hours"))
  )

########## Final Table for Publication ##############
Table_GC <- GC_filled %>%
  select(Total_Time,Amount,Setup,Replicate) %>%
  as.data.frame()
#####################################################


# # Fiding the Best Model Fit per Setup and Replicate
# fit_models_per_group <- Table_GC %>%
#   group_by(Setup, Replicate) %>%
#   group_modify(~ {
#     df <- .x
#     
#     # Prepare safe fitting functions that return NA if they fail
#     fit_linear <- tryCatch(
#       lm(Amount ~ Total_Time, data = df),
#       error = function(e) NA
#     )
#     
#     fit_exp <- tryCatch({
#       # Only fit if all Amount > 0 (log needed)
#       if (all(df$Amount > 0)) {
#         lm(log(Amount) ~ Total_Time, data = df)
#       } else {
#         NA
#       }
#     }, error = function(e) NA)
#     
#     fit_sigmoid <- tryCatch({
#       nls(
#         Amount ~ a / (1 + exp(-(Total_Time - x0) / b)),
#         data = df,
#         start = list(
#           a = max(df$Amount, na.rm = TRUE),
#           x0 = mean(df$Total_Time, na.rm = TRUE),
#           b = 1
#         ),
#         control = nls.control(warnOnly = TRUE)
#       )
#     }, error = function(e) NA)
#     
#     # Return models as a list-column tibble
#     tibble(
#       model_linear = list(fit_linear),
#       model_exp = list(fit_exp),
#       model_sigmoid = list(fit_sigmoid)
#     )
#   }) %>%
#   ungroup()
# 
# best_model_stats <- fit_models_per_group %>%
#   mutate(
#     # Get AIC or NA if model failed
#     aic_linear = map_dbl(model_linear, ~ if(inherits(.x, "lm")) AIC(.x) else NA_real_) ,
#     aic_exp    = map_dbl(model_exp,    ~ if(inherits(.x, "lm")) AIC(.x) else NA_real_) ,
#     aic_sigmoid= map_dbl(model_sigmoid,~ if(inherits(.x, "nls")) AIC(.x) else NA_real_)
#   ) %>%
#   rowwise() %>%
#   mutate(
#     # Find model name with lowest AIC
#     best_model = c("linear", "exp", "sigmoid")[which.min(c(aic_linear, aic_exp, aic_sigmoid))]
#   ) %>%
#   ungroup() %>%
#   select(Setup, Replicate, aic_linear, aic_exp, aic_sigmoid, best_model)

# Best fit model search for each setup
# 1. Filter data (remove 'B' replicates and 'Blank' Setup)
filtered_data <- Table_GC %>%
  filter(!is.na(Replicate)) %>%
  filter(Replicate != "B", Setup != "Blank") #%>%
  # Manaully Removed due to potentially logarithmic
  filter(!(Replicate == "3" & Setup == "Control" & Total_Time == 73.5) ) %>%
  filter(!(Replicate == "3" & Setup == "Control" & Total_Time == 97) ) 
  
  
# 2. Define fitting function (same as before)
fit_models <- function(df) {
  start_linear <- list(a = min(df$Amount), b = 0)
  start_exp <- list(a = min(df$Amount)+1e-6, b = 0.01)
  start_sigmoid <- list(a = max(df$Amount), b = 0.1, c = median(df$Total_Time), d = min(df$Amount))
  
  fit_linear <- try(lm(Amount ~ Total_Time, data = df), silent = TRUE)
  fit_exp <- try(nlsLM(Amount ~ a * exp(b * Total_Time),
                       data = df,
                       start = start_exp,
                       control = nls.lm.control(maxiter=100)),
                 silent = TRUE)
  fit_sigmoid <- try(nlsLM(Amount ~ a / (1 + exp(-b*(Total_Time - c))) + d,
                           data = df,
                           start = start_sigmoid,
                           control = nls.lm.control(maxiter=200)),
                     silent = TRUE)
  
  aic_linear <- if(inherits(fit_linear, "lm")) AIC(fit_linear) else NA
  aic_exp <- if(inherits(fit_exp, "nls")) AIC(fit_exp) else NA
  aic_sigmoid <- if(inherits(fit_sigmoid, "nls")) AIC(fit_sigmoid) else NA
  
  aic_vals <- c(linear = aic_linear, exponential = aic_exp, sigmoid = aic_sigmoid)
  best_model <- names(which.min(aic_vals))
  
  list(
    fit_linear = fit_linear,
    fit_exp = fit_exp,
    fit_sigmoid = fit_sigmoid,
    best_model = best_model,
    aic = aic_vals
  )
}

# 2.5 Setting up Corrected Data 
useful_data <-  filtered_data #%>%
  # Manually set the linear fitting for those above this value
  #filter(Setup == 'Coal' | (Setup != "Coal" & Amount > 1 ) ) %>%  For all Days with Coal as no growth
  #filter( !(Total_Time >= 97  & Setup %in% c('Shale','Sand')) ) # For Day 36
  filter( !(Total_Time > 100 & Setup %in% c('Shale','Sand','Coal')) ) %>% filter(Amount > 1) # For Day 29 
  

# 3. Fit models by Setup and keep the best fit object
fit_results <- useful_data %>%
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
                            linear = predict(fits$fit_linear, newdata = data.frame(Total_Time = new_time)),
                            exponential = {
                              a <- coef(fits$fit_exp)["a"]
                              b <- coef(fits$fit_exp)["b"]
                              a * exp(b * new_time)
                            },
                            sigmoid = {
                              coefs <- coef(fits$fit_sigmoid)
                              a <- coefs["a"]
                              b <- coefs["b"]
                              c <- coefs["c"]
                              d <- coefs["d"]
                              a / (1 + exp(-b * (new_time - c))) + d
                            })
      
      tibble(
        Total_Time = new_time,
        Amount = pred_amount
      )
    })
  ) %>%
  select(Setup, pred_df) %>%
  unnest(pred_df)

# 5. Plot with ggplot2 - data points + best fit lines, facet by Setup
#################################### PLot for Pub ############################################
Plot_DNP <- ggplot() +
  geom_point(data = filtered_data, aes(x = Total_Time, y = Amount), color = "Black") +
  geom_line(data = predictions, aes(x = Total_Time, y = Amount), color = "red", size = 1) +
  facet_wrap(~Setup) +
  labs(title = paste0("Day ",input_day),
       x = "Total Time",
       y = "Amount") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.placement = "outside",
    strip.background = element_rect(fill="grey90"),
    strip.text = element_text(face = "bold", size = 18),
    #strip.switch.pad.wrap = unit(0.2, "cm"),
    legend.position = "none"
  )

Plot_DNP
#################################### PLot for Pub ############################################

#6. Compute the Slopes and save into a table
slopes_table <- useful_data %>%
    group_by(Setup, Replicate) %>%
  summarise(
    # Fit a linear model Amount ~ Total_Time for each group
    slope = if(n() > 1) {
      coef(lm(Amount ~ Total_Time))[2]  # slope coefficient
    } else {
      NA_real_  # Not enough data points to fit
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
    Replicates = n()
            ) %>%
  mutate(DAY = input_day) # MANUALLY SET THE DAY in line
#################################################


# 7. Alternative output format for plot
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

# Create the table grob with the minimal theme
inset_table <- Table_Slopes %>%
  unite('Slope', Mean, Std_Dev, sep = '±') %>%
  select(Setup, Slope)
table_grob <- tableGrob(
  inset_table, 
  rows = NULL, theme = minimal_theme
  )
Plot_DNP_Table <- Plot_DNP + inset_element(
  table_grob,
  left = 0.75, bottom = 0.10, right = .95, top = .35
)
Plot_DNP_Table



################ SAVE the Tables and Plots for Publication#####################################
ggsave(plot = Plot_DNP, filename = paste0("~/Workbench/P3/Expt2/Plot_DNP_",input_day,".png"),
       units = c('px'),
       width = 2000, height = 1200, dpi = 120)

ggsave(plot = Plot_DNP_Table, filename = paste0("~/Workbench/P3/Expt2/Plot_DNP_InsetTable",input_day,".png"),
       units = c('px'),
       width = 2000, height = 1200, dpi = 120)

write_csv(Table_GC, file = paste0("~/Workbench/P3/Expt2/Table_GC_Summary_",input_day,".csv") )

write_csv(Table_Slopes, file = paste0("~/Workbench/P3/Expt2/Table_GC_Slopes_",input_day,".csv") )
################ SAVE the Tables and Plots for Publication#####################################

