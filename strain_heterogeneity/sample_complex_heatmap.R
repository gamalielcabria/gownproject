# Initialize Libraries
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

# Read the FastANI matrix file using read.csv()
input_file<-"/home/glbcabria/Workbench/Aquificales/02-ANI_CA_NZ/ALLNZ_UBA11096_75UP_CA/fastANI.out.matrix"

# Read the first line to get the number of genomes
num_genomes <- as.numeric(readLines(input_file, n = 1))

# Read the rest of the file into a single column
ani_data <- read_delim(input_file, delim = ';', col_names = FALSE, skip = 1 )

# Extract genome names (first column)
genome_names <- gsub(".fasta.*$", '', ani_data[[1]])

# Convert the column into a dataframe with edited column and row names
ani_matrix <- ani_data %>%
  separate(X1, into=c('row_id',genome_names), sep = '\t') %>%
  mutate(row_id = gsub('.fasta','',row_id)) %>%
  column_to_rownames('row_id') %>%
  mutate(across(everything(),as.numeric)) %>%
  as.matrix()

# Fill diagonal with 100% ANI (self-comparison)
diag(ani_matrix) <- 100

# Mirror the lower triangle to the upper triangle
ani_matrix[upper.tri(ani_matrix)] <- t(ani_matrix)[upper.tri(ani_matrix)]




# Plot the data

# Define color function based on the given ranges
col_fun <- colorRamp2(
  c(70.00, 86.00, 95.99, 96.00, 99.29, 99.30, 99.80, 100.00), 
  c("#4575b4", "#66ad75", "#b1d081", "#ffffff","#fff197","#f9b968", "#ed7e51", "#d43d51")
)

# Create the heatmap
plot <- Heatmap(ani_matrix,
        col = col_fun, 
        name = "Value", 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_names = TRUE, 
        show_column_names = TRUE,
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 14),  # Increase legend title font size
          labels_gp = gpar(fontsize = 12), # Increase legend labels font size
          legend_height = unit(6, "cm"),   # Increase legend height
          legend_width = unit(1, "cm"),    # Increase legend width
          at = c(70,85, 90, 96, 99.3, 99.8, 100),       # Add demarcation at 85, 96, and 99.3
          labels = c("70","85", "90","96", "99.3", "99.8", "100")  # Explicit labels for demarcations
        )
        )
plot
