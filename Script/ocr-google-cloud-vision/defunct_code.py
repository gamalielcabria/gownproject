### GOOGLE OCR ###

from google.cloud import vision

# Initialize Google Vision API client

def perform_ocr(image_path):
    vision_client = vision.ImageAnnotatorClient()
    with open(image_path, "rb") as image_file:
        content = image_file.read()
    image = vision.Image(content=content)
    response = vision_client.document_text_detection(image=image)
    full_text_annotation = response.full_text_annotation

    if full_text_annotation:
        lines = []
        for page in full_text_annotation.pages:
            for block in page.blocks:
                for paragraph in block.paragraphs:
                    line_text = ""
                    for word in paragraph.words:
                        word_text = "".join([symbol.text for symbol in word.symbols])
                        line_text += word_text + " "
                    lines.append(line_text.strip())
        return lines  # Return each line as a separate string in a list
    
    if response.error.message:
        raise Exception(f"Vision API error: {response.error.message}")
    return []

# # Perform OCR
# extracted_lines = perform_ocr('compressed_image.jpg')
# if not extracted_lines:
#     print(f"Warning: No text detected in {image_file}.")
#     continue
# print(f"File: {image_path} \n Text: {extracted_lines}")
# # Perform OCR
# extracted_lines = perform_ocr('compressed_image.jpg')
# if not extracted_lines:
#     print(f"Warning: No text detected in {image_file}.")
#     continue
# print(f"File: {image_path} \n Text: {extracted_lines}")

### GOOGLE OCR ###