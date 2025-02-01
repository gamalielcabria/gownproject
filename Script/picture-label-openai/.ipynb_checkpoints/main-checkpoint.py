import os, re
from openai import OpenAI
from utility_code import encode_image
from utility_code import extract_details
from PIL import Image


client = OpenAI()

folder_path = '/home/glbcabria/Workbench/GOWN_Data/Photos/GOWN_Bottles'

# if os.path.isdir(input_dir):
#     folder_path=input_dir
# else:
#     folder_path=os.getcwd

for filenames in os.listdir(folder_path):
    if filenames.endswith(',jpeg') or filenames.endswith('.jpg'):
        print(filenames)
        image_file = filenames
        image_path = os.path.join(folder_path,image_file)

        ## Resizing input image
        image = Image.open(image_path)
        width, height = image.size
        new_size = (width//1, height//1)
        resized_image = image.resize(new_size)
        resized_image.save('compressed_image.jpg', optimize=True, quality=50)

        base64_image = encode_image('compressed_image.jpg')
        # print(f"{base64_image}\n")
        
        response = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[
                {
                    "role": "user",
                    "content": [
                        {
                            "type": "text",
                            "text": "The label in the images can be curved and the text therefor can appear contorted or blurred. Can you convert the label and the handwriting in this picture to text? Please identify date and time separately. Also, identify texts with regex pattern '(2\d+GW.*)' as sample number. Also, identify location, site, or station if available. Please output date as YYYY-MM-DD. Also, identify a standalone number that is part of the location or site label as 'ID', if available. The ID can be a standalone 3 or 4 digit number. Please output a text in the format of a key value pair."
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


        ## Responses:
        # print('Completion tokens: ', response.usage.completion_tokens)
        # print('Prompt tokens: ', response.usage.prompt_tokens)
        # print('Total tokens: ', response.usage.total_tokens)
        # print(response.choices[0].message)
        # print(response.choices[0].message.content)

        details = extract_details(response.choices[0].message.content)
        #print(details)

        pattern_two = r'(Sample Number:.*)'
        if details:
            sample_number_match = re.search(pattern_two, details)
            if sample_number_match:
                sample_number = sample_number_match.group().split(': ')[1]
            else:
                sample_number = f'No_Sample_Number'
            print(sample_number)

            if '/' in sample_number:
                new_sample_number = sample_number.replace('/','-')
                new_txt_filename = f'{new_sample_number}.gpt4.txt'
                new_image_filename = f'{new_sample_number}_gpt4.jpg'
            else:
                new_txt_filename = f'{sample_number}.gpt4.txt'
                new_image_filename = f'{sample_number}_gpt4.jpg'

            with open(new_txt_filename, "w") as file:
                file.write(details)
            
            resized_image.save(new_image_filename, optimize=True, quality=50)
        else:
            print(f'No Details Exist for {image_file}')