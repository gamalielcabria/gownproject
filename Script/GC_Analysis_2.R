#install.packages("pdftools")
library(pdftools)
library(tidyverse)
library(broom)
setwd('/home/glbcabria/Workbench/P3/Expt2/GC')

##################### Functions #############################################
find_N2O_firstmatch <- function(text, lower = 1.8, upper = 1.9) {
  # Split on literal "\n"
  txt_lines <- strsplit(text, "\\\n")[[1]]
  
  all_numbers <- regmatches(txt_lines, gregexpr("\\d+\\.\\d+", txt_lines))
  
  for (i in seq_along(all_numbers)) {
    nums <- as.numeric(all_numbers[[i]])
    match_nums <- nums[nums >= lower & nums <= upper]
    if (length(match_nums) > 0) {
      return(line = txt_lines[i])  #number = match_nums[1])
    }
  }
  return(NULL)
}
#############################################################################

# File Location
files <- list.files(pattern="\\.pdf$", ignore.case = TRUE)


# Create empty dataframe to store results
GC_df <- data.frame(File = character(),
                         Results = character(),
                         stringsAsFactors = FALSE)

# Read through files
for (file in files) {
  txt <- pdf_text(file)
  #filename <- regmatches(txt[2], regexpr("\\d{8}[-_][A-Za-z0-9_-]+\\.D", txt[2]))
  filename <- gsub("\\.pdf",'',as.character(file))
  N2O_line <- trimws(find_N2O_firstmatch(txt[2]))
  
  results <- gsub("\\s+", "_", N2O_line)
  
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
                       as.numeric(map_chr(str_split(Results, "_"), ~ .x[length(.x) -1]) ),
                       as.numeric(map_chr(str_split(Results, "_"), ~ .x[length(.x) -3]) ) 
                       ) 
         )%>%
  filter(!grepl("PREPROCESS",File, ignore.case = TRUE)) %>%
  mutate(File = gsub("-POSTPROCESS", "", File)) %>%
  separate(File, into = c('Date-Time','Type','Sample'), sep = "_") %>%
  select(-Results) %>%
  separate('Date-Time', into = c("Date","Time"), sep = '-' ) %>%
  #mutate(Amount = ifelse(Type != 'Calibration', NA, Amount) ) 
  mutate(Amount = ifelse(grepl('Standards-0', Sample), 0, Amount ),
         Area = ifelse(grepl('Standards-0', Sample), 0, Area ) ) %>%
  mutate(Sample = gsub("-\\d+$|-\\d{1}[A-Za-z]{2}$", "", Sample) ) %>% # Remove the -8PM/-1/-2 from Sample Name
  replace_na(list(Area = 0))

# Calibration Model
models <- GC_processed_df %>%
  filter(Type == 'Calibration') %>%
  group_by(Date) %>%
  summarise(
    model = list(lm(Amount ~ Area, data = cur_data()))
  )

# Join models back to main dataframe by Date
GC_filled <- GC_processed_df %>%
  left_join(models, by = "Date") %>%
  rowwise() %>%
  mutate(
    # If Type is Sample and Amount is NA, predict Amount using the model and Area
    Amount = if_else(
      Type == "Sample" & is.na(Amount),
      predict(model[[1]], newdata = tibble(Area = Area)),
      Amount
    )
  ) %>%
  ungroup() %>%
  select(-model)  # remove model column
