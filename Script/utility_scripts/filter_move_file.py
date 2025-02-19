import pandas as pd
import shutil
import os

# Define input, output file paths, and destination folder
input_file = "~/Workbench/Aquificales/Genomes/CA_genomes/checkm2_out/quality_report.tsv"  # Change this to the actual file path
output_file = "filtered_output.txt"
destination_folder = "/home/glbcabria/Workbench/Aquificales/01-PopCOGenT_CA_NZ/"  # Change this to the desired folder

# Ensure the destination folder exists
os.makedirs(destination_folder, exist_ok=True)

# Load the tab-delimited file into a pandas DataFrame
df = pd.read_csv(input_file, sep='\t')

# Filter rows where 'Completeness' is greater than 50
filtered_df = df[df['Completeness'] > 75]

# Save the filtered results to a new file
filtered_df.to_csv(output_file, sep='\t', index=False)

# Move selected files to the destination folder
for name in filtered_df['Name']:
    source_path = f"{name}.fna"  # Assuming files are named based on the 'Name' column and the script is run where the files are
    destination_path = os.path.join(destination_folder, f"{name}_CA.fasta")
    if os.path.exists(source_path):
        shutil.copy(source_path, destination_path)

print(f"Filtered data saved to {output_file}")
print(f"Selected files moved to {destination_folder}")