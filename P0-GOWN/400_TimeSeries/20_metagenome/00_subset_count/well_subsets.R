# =========================================================
# Count all shared year subsets across wells
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
})

# ---------------------------------------------------------
# 1) Input
# ---------------------------------------------------------
setwd("/home/gam/github/gownproject/P0-GOWN/300_TimeSeries/00_subset_count")
input_file <- "/home/gam/github/gownproject/P0-GOWN/300_TimeSeries/04_well_presence_absence.csv"

df <- read_csv(input_file, show_col_types = FALSE)

# ---------------------------------------------------------
# 2) Checks and cleanup
# ---------------------------------------------------------
if (!"Well" %in% names(df)) {
  stop("Input must contain a column named 'Well'")
}

year_cols <- setdiff(names(df), "Well")

if (length(year_cols) < 3) {
  stop("Need at least 3 year columns to make 3-year subsets")
}

df <- df %>%
  mutate(across(all_of(year_cols), as.numeric)) %>%
  mutate(across(all_of(year_cols), ~ tidyr::replace_na(., 0)))

bad_vals <- unlist(df[year_cols], use.names = FALSE)
if (any(!bad_vals %in% c(0, 1))) {
  stop("Year columns must contain only 0/1/NA")
}

# ---------------------------------------------------------
# 3) Helper: generate all subsets for one well
# ---------------------------------------------------------
# For a well present in years c("2016","2019","2022","2024"),
# this returns all combinations of size 3 up to all years present.
get_all_subsets <- function(well_name, present_years, min_size = 3) {
  
  n <- length(present_years)
  
  if (n < min_size) {
    return(NULL)
  }
  
  subset_list <- lapply(min_size:n, function(k) {
    combos <- combn(present_years, k, simplify = FALSE)
    
    tibble(
      Well = well_name,
      subset_size = k,
      subset = vapply(combos, paste, collapse = "_", FUN.VALUE = character(1))
    )
  })
  
  bind_rows(subset_list)
}

# ---------------------------------------------------------
# 4) Expand every well into all possible subsets
# ---------------------------------------------------------
all_subsets_long <- map_dfr(seq_len(nrow(df)), function(i) {
  
  row_i <- df[i, ]
  well_i <- row_i$Well
  
  present_years <- year_cols[as.numeric(row_i[1, year_cols]) == 1]
  
  get_all_subsets(
    well_name = well_i,
    present_years = present_years,
    min_size = 2
  )
})

# ---------------------------------------------------------
# 5) Count how many wells support each subset
# ---------------------------------------------------------
subset_counts <- all_subsets_long %>%
  distinct(Well, subset_size, subset) %>%
  count(subset_size, subset, name = "n_wells") %>%
  arrange(desc(subset_size), desc(n_wells), subset)

# ---------------------------------------------------------
# 6) List which wells support each subset
# ---------------------------------------------------------
subset_members <- all_subsets_long %>%
  distinct(Well, subset_size, subset) %>%
  arrange(desc(subset_size), subset, Well) %>%
  group_by(subset_size, subset) %>%
  summarise(
    n_wells = n(),
    wells = paste(sort(Well), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(subset_size), desc(n_wells), subset)

# ---------------------------------------------------------
# 7) Summary by subset size only
# ---------------------------------------------------------
# Example:
# among all 3-year subsets, how many unique subsets exist,
# and what is the max support?
subset_size_summary <- subset_counts %>%
  group_by(subset_size) %>%
  summarise(
    n_unique_subsets = n(),
    max_wells_in_a_subset = max(n_wells),
    .groups = "drop"
  ) %>%
  arrange(desc(subset_size))

# ---------------------------------------------------------
# 8) Top subsets within each size
# ---------------------------------------------------------
top_subsets_by_size <- subset_members %>%
  group_by(subset_size) %>%
  filter(n_wells == max(n_wells)) %>%
  ungroup() %>%
  arrange(desc(subset_size), subset)

# ---------------------------------------------------------
# 9) Optional: count number of years per well too
# ---------------------------------------------------------
well_year_counts <- df %>%
  mutate(n_years = rowSums(across(all_of(year_cols)))) %>%
  arrange(desc(n_years), Well)

# ---------------------------------------------------------
# 10) Print results
# ---------------------------------------------------------
cat("\n====================================================\n")
cat("Top of all subset counts\n")
cat("====================================================\n")
print(subset_counts, n = 50)

cat("\n====================================================\n")
cat("Subset members\n")
cat("====================================================\n")
print(subset_members, n = 50, width = Inf)

cat("\n====================================================\n")
cat("Subset-size summary\n")
cat("====================================================\n")
print(subset_size_summary, n = Inf)

cat("\n====================================================\n")
cat("Most supported subset(s) within each subset size\n")
cat("====================================================\n")
print(top_subsets_by_size, n = Inf, width = Inf)

# ---------------------------------------------------------
# 11) Export results
# ---------------------------------------------------------
write_tsv(all_subsets_long,      "01_all_well_supported_subsets_long.tsv")
write_tsv(subset_counts,         "02_all_subset_counts.tsv")
write_tsv(subset_members,        "03_all_subset_members.tsv")
write_tsv(subset_size_summary,   "04_subset_size_summary.tsv")
write_tsv(top_subsets_by_size,   "05_top_subsets_by_size.tsv")
write_tsv(well_year_counts,      "06_well_year_counts.tsv")

cat("\nDone. Files written:\n")
cat("  01_all_well_supported_subsets_long.tsv\n")
cat("  02_all_subset_counts.tsv\n")
cat("  03_all_subset_members.tsv\n")
cat("  04_subset_size_summary.tsv\n")
cat("  05_top_subsets_by_size.tsv\n")
cat("  06_well_year_counts.tsv\n")