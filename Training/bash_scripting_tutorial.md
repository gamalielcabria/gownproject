
---
title: Bash Scripting Tutorial
nav_order: 4
parent: Intro to Bioinfo
---

# Bash Scripting Tutorial

This tutorial introduces the basics of Bash scripting, useful for automating tasks in bioinformatics.

---

## 🧩 Basic Bash Scripting

A Bash script is a plain text file with commands executed sequentially by the Bash shell.

### 📝 Creating a Script

```bash
nano my_script.sh
```

Add the following:

```bash
#!/bin/bash
echo "Hello from your first script!"
```

Then make it executable and run it:

```bash
chmod +x my_script.sh
./my_script.sh
```

---

## 🛠️ Variables and Input

Variables help reuse values and pass arguments.

```bash
#!/bin/bash
sample="Sample1"
echo "Processing data for $sample"
```

**Using arguments**:

```bash
#!/bin/bash
echo "This script was called with: $1 and $2"
```

Run with:

```bash
./my_script.sh file1.txt file2.txt
```

---

## 🔁 Loops

### For Loops

```bash
for file in *.fastq; do
    echo "Processing $file"
done
```

### While Loops

```bash
count=1
while [ $count -le 5 ]; do
  echo "Count is $count"
  ((count++))
done
```

---

## ⚙️ Conditional Statements

```bash
if [ -f "genome.fasta" ]; then
    echo "Genome file found!"
else
    echo "Genome file is missing!"
fi
```

---

## Real-World Example

```bash
#!/bin/bash

indir="raw_data"
outdir="qc_results"
mkdir -p "$outdir"

for file in "$indir"/*.fastq; do
    if [ -s "$file" ]; then
        echo "Running FastQC on $file"
        fastqc "$file" -o "$outdir"
    else
        echo "Skipping $file, it is empty."
    fi
done

echo "All files processed."
```

---

## Try It Yourself!

Create a script that checks if a `.gz` file exists, unzips it, and prints how many lines it has.

```bash
#!/bin/bash
file="example.gz"
if [ -f "$file" ]; then
    echo "Unzipping $file"
    gunzip -c "$file" | wc -l
else
    echo "$file not found."
fi
```

---

## 🔚 Summary

You’ve now learned how to:
- Create and execute bash scripts
- Use variables and arguments
- Run loops and conditionals
- Process real data automatically

Enjoy scripting!

