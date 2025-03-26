from utility_code import encode_image, extract_details
import os
import pandas as pd
from PIL import Image
from openai import OpenAI
from datetime import datetime
 

def gc_img_openai_processing(folder_path, output_file=None):
    # Initialize OpenAI API client
    openai_client = OpenAI()

    # Validate folder_path
    if not os.path.exists(folder_path) or not os.path.isdir(folder_path):
        print(f"Invalid folder path: {folder_path}")
        return

    # Data storage for pandas later
    data = []


    # Process each image in the folder
    for filenames in os.listdir(folder_path):
        if filenames.endswith('.jpeg') or filenames.endswith('.jpg'):
            print(f"Processing file: {filenames}")
            image_file = filenames
            image_path = os.path.join(folder_path, image_file)

            # Resize image
            image = Image.open(image_path)
            width, height = image.size
            resized_image = image.resize((width, height))
            resized_image.save('compressed_image.jpg', optimize=True, quality=50)
            base64_image = encode_image('compressed_image.jpg')
            

            # Send extracted text to OpenAI API
            try:
                response = openai_client.chat.completions.create(
                    model="gpt-4o-mini",
                    messages=[
                        {
                            "role": "user",
                            "content": [
                                {
                                    "type": "text",
                                    "text": "Extract the table the text from the image. The table has column names: RetTime [min], Type, Area [Hz*s], Amt/Area, Amount [Mole %], Grp, Name . The rows has values like N2O and H2S. Add a column with values of the filename only found in the end of the line starting from Data File or found with the extension .D, do not include the .D extension in the output. Do not include to the table the column for Name, Type, and Grp. Output it as a code block of csv formatted table."
                                },
                                {
                                    "type": "image_url",
                                    "image_url": {
                                        "url": f'data:image/jpeg;base64,{base64_image}',
                                        "detail": "high"
                                    }
                                }
                            ]
                        }
                    ] #,max_tokens=300
                )
                # Parse OpenAI API response
                # print(response.choices[0].message)
                details = extract_details(response.choices[0].message.content)
                #print(f"extracted details: {details}")
                line = details.strip().split('\n')[2]
                array = line.split(",")
                array.append(os.path.splitext(filenames)[0])
                if len(array) == 6:
                    data.append(array)
                else:
                    print(f"Warning: File {image_file} does not contain exactly 5 elements in the first line.")
            
            except Exception as e:
                print(f"Error processing {image_file}: {e}")

    # Save results to CSV
    if data:
        df = pd.DataFrame(data, columns=['RetTime','Area','Amt/Area','Amount','Filename','ImageName'])
        current_date = datetime.now().strftime("%Y-%m-%d_%H%M")
        
        if output_file is None:
            output_file = f"{folder_path}/GC_output_gpt_{current_date}.csv"
        try:
            df.to_csv(output_file, index=False)
            print(f"Results saved to {output_file}")
        except Exception as e:
            print(f"Error saving file: {e}")
        return df
    else:
        print("No valid data extracted")


# # # Path to the folder containing images for testing purposes
folder_path = '/home/glbcabria/Workbench/P3/Results/GCPhotos/GWMix'
gc_img_openai_processing(folder_path)


