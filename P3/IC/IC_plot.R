#R
library(tidyverse)

folder <- "/home/gam/github/gownproject/P3/IC"
setwd(folder)
files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)

for (file in files) {
  # Create a safe variable name based on file name
  var_name <- tools::file_path_sans_ext(basename(file))
  var_name <- make.names(paste0("DF_",var_name))  # make it a valid R variable name
  
  # Read and clean
  lines <- readLines(file)
  lines <- lines[-c(1:10, 12, 13)]
  
  df <- read.csv(text = lines, sep = ",", header = TRUE, 
                   fill = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  # Drop "No." column if present
  if ("No." %in% names(df)) df <- df[, names(df) != "No."]

  # Assign to a variable with the file name
  area_row2 <- if ("Area" %in% names(df))
                 suppressWarnings(df$Area[1]) else NA_real_
  df$Chemical <- rep(area_row2, nrow(df))
  df <- df[-1,]
  df <- df %>%
    mutate(Injection.Name = gsub('_PT','', Injection.Name)) %>%
    filter(!(grepl("STD",Injection.Name))) %>%
    mutate(Timepoint = str_extract(Injection.Name, regex("(?<![A-Za-z0-9])D\\d+(?![A-Za-z0-9])", ignore_case = TRUE) )) %>%
    mutate(Sample = str_extract(Injection.Name, "(?<=_)[[:alpha:]]+(?=\\d*$)")) %>%
    filter(Sample != "<NA>") %>%  
    mutate(Amount = ifelse(Amount == "n.a.", 0, Amount)) %>%
    mutate(Amount = as.numeric(Amount)) %>%
    select(Sample, Timepoint, Amount, Chemical)


  assign(var_name, df)
  cat("Loaded:", var_name, " | Area_row2 =", area_row2, "\n")
}

# collect all objects whose names start with DF_
df_names <- ls(pattern = "^DF_")          # (fixed from "^DF_$")
if (!length(df_names)) stop("No DF_* data frames found.")

dfs <- mget(df_names, inherits = TRUE)
dfs <- dfs[vapply(dfs, is.data.frame, logical(1))]

# combine; keeps existing Source column from each DF_*
dat <- bind_rows(dfs)

# Plot the values
dat_num <- dat %>%
  mutate(
    Timepoint = readr::parse_number(Timepoint),
    Amount    = as.numeric(Amount),
    Sample    = recode(Sample,
      CO = "Coal",
      CL = "Control",
      SH = "Shale",
      SS = "Sandstone",
      SA = "Sand",
      MILLIQ = "BLANK",
      ULTRA = "BLANK",
      DNU = "BLANK"
    )
  )

sum_df <- dat_num %>%
  filter(Sample != "BLANK") %>%
  group_by(Chemical, Sample, Timepoint) %>%
  summarise(
    n    = n(),
    mean = mean(Amount, na.rm = TRUE),
    sd   = sd(Amount,   na.rm = TRUE),
    .groups = "drop"
  )


library(ggsci)


plot_IC <- ggplot(sum_df, aes(x = Timepoint, y = mean, color = Sample, group = Sample)) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 3.6) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.15, linewidth = 0.8) +
  facet_wrap(~ Chemical, scales = "free_y") +
  scale_color_jco(name = "Sample") +         # use jco palette for lines/legend
  labs(x = "Timepoint (days)", y = "Mean Amount (mg/L)") +
  theme_minimal(base_size = 13)+
  theme(
    strip.text = element_text(size = 10, face = "bold"),   # small + minimal
    strip.background = element_blank()                     # removes shading/box
  )

ggsave("IC_plot_day78.png", plot_IC, width = 1800, height = 1000, units = "px", dpi = 150)
#write.csv(sum_df, "IC_summary_leaching_data.csv", row.names = FALSE)
#write.csv(dat_num, "IC_all_leaching_data.csv", row.names = FALSE)
