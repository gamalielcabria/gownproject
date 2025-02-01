import pandas as pd

# Example extracted text (from Vision API or OCR)
extracted_text = [
    '1 1/11/2025 6:06:12 PM SYSTEM',
    'Data File C : \\ CHEM32 \\ 1 \\ DATA \\ GAMA \\ N20 20250111_551PM 100N200ARCAL CURVE.D Sample Name : D2 headspace',
    'RetTime Type',
    '[ min ]',
    'Area [ Hz * s ]',
    'Amt / Area',
    '- 1',
    '-- | ---',
    '---------',
    '1.323 BB S',
    'Amount [ Mole % ] | ---------- | -- | -- 6.75129e4 1.52111e - 6 1.02695e - 1',
    'Grp',
    'Name',
    'N20',
    '2.091',
    'H2S',
    'Totals :',
    'Signal 2 : FID1 A , Front Signal',
    'SHOT ON POCO ReNFame Type',
    'Area',
    '[ min ]',
    '[ pA * s ]',
    '1.02695e - 1',
    'Amt / Area',
    'Amount',
    'Grp Name',
    '[ Mole % ]',
    '--------',
    'Page',
    '1 of 2'
]

# Process the extracted text to extract relevant table rows
def parse_text_to_table(extracted_text):
    table_data = []
    for i, line in enumerate(extracted_text):
        if '1.323' in line or '2.091' in line:  # Detect table rows
            # Clean and split row content
            parts = line.split()
            if len(parts) >= 7:
                table_data.append({
                    "RetTime [min]": parts[0],
                    "Type": parts[1],
                    "Area [Hz*s]": parts[2],
                    "Amt/Area": parts[3],
                    "Amount [Mole %]": parts[4],
                    "Grp": parts[5],
                    "Name": parts[6]
                })
    return table_data

# Extract table data
table_data = parse_text_to_table(extracted_text)

# Create a DataFrame
df = pd.DataFrame(table_data)

# Display the DataFrame
print(df)

# Save to a CSV file
#df.to_csv("extracted_table.csv", index=False)
