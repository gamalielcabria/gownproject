import base64
import re

"""
Utility code fetch from Jie Jenn's Image analysis With OpenAI API GPT4 Vision In Python
(https://www.youtube.com/watch?v=_1ujhANv6a4)
"""


def encode_image(image_path):
    """
    Encodes and image to base 64 String
    """
    with open(image_path, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode('utf-8')

def extract_code(text):
    # Pattern to match the code block
    pattern = r'```python\s+(.*?)\s+```'

    #extract the code
    match = re.search(pattern, text, re.DOTALL)
    if match:
        return match.group(1)
    return None

def extract_details(text):
    pattern_one = r'```(.*?)```'
    match = re.search(pattern_one,text, re.DOTALL)
    if match:
        details = match.group(1)
        return details
    return None