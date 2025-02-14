import os
import pillow_heif
from PIL import Image
import concurrent.futures

# Function to convert a single HEIC file to JPEG
def convert_heic_to_jpeg(heic_path, output_dir):
    try:
        # Open HEIC file
        heif_file = pillow_heif.open_heif(heic_path)

        # Extract image data properly
        image = Image.frombytes(
            heif_file.mode, 
            heif_file.size, 
            heif_file.data.tobytes()  # Ensures correct byte format
        )

        # Define output filename
        output_filename = os.path.splitext(os.path.basename(heic_path))[0] + ".jpg"
        output_path = os.path.join(output_dir, output_filename)

        # Save as JPEG with quality settings
        image.save(output_path, "JPEG", quality=90)
        print(f"Converted: {heic_path} -> {output_path}")
    
    except Exception as e:
        print(f"Error processing {heic_path}: {e}")

# Function to batch convert all HEIC files in a folder using multi-threading
def batch_convert_heic_to_jpeg(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    heic_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.lower().endswith(".heic")]

    # Use multi-threading to speed up conversion
    with concurrent.futures.ThreadPoolExecutor() as executor:
        executor.map(lambda f: convert_heic_to_jpeg(f, output_dir), heic_files)

    print("Batch conversion complete!")

# Example usage
input_directory = "/home/glbcabria/Workbench/P3/Results/GCPhotos/D2_Madison"  # Change this to your HEIC folder
output_directory = "/home/glbcabria/Workbench/P3/Results/GCPhotos/D2_Madison/JPEGs"  # Change to your output folder

batch_convert_heic_to_jpeg(input_directory, output_directory)
