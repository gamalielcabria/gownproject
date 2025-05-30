
# 🧪 Bioinformatics Practice Exercise: Linux Command Line

This exercise uses files from [this GitHub folder](https://github.com/gamalielcabria/gamalielcabria.github.io/tree/main/Intro2Bioinfo/files).

To get started, download the files to your terminal using `wget`:

```bash
mkdir -p ~/bioinfo_practice && cd ~/bioinfo_practice

wget -P ./ https://raw.githubusercontent.com/gamalielcabria/gamalielcabria.github.io/main/Intro2Bioinfo/files/data.tsv
wget -P ./ https://raw.githubusercontent.com/gamalielcabria/gamalielcabria.github.io/main/Intro2Bioinfo/files/gene_counts.tsv
wget -P ./ https://raw.githubusercontent.com/gamalielcabria/gamalielcabria.github.io/main/Intro2Bioinfo/files/genes.txt
wget -P ./ https://raw.githubusercontent.com/gamalielcabria/gamalielcabria.github.io/main/Intro2Bioinfo/files/lengths.txt
wget -P ./ https://raw.githubusercontent.com/gamalielcabria/gamalielcabria.github.io/main/Intro2Bioinfo/files/random.fasta
```

---

## 📝 Exercises

### 1. **Inspecting Files**
- View the first few lines of each file using `head`.
- Count the number of rows and columns in `gene_counts.tsv`.

```bash
head gene_counts.tsv
wc -l gene_counts.tsv
head -1 gene_counts.tsv | awk -F'\t' '{print NF}'
```

---

### 2. **Extracting and Combining Columns**
- Extract the 1st and 3rd columns from `data.tsv` and save to `score.tsv`.
- Extract the 1st and 2nd columns and save to `value.tsv`.
- Merge them side by side using `paste`.

```bash
cut -f1,3 data.tsv > score.tsv
cut -f1,2 data.tsv > value.tsv
paste score.tsv value.tsv > combined.tsv
```

---

### 3. **Grep and Regex**
- Count how many genes start with `geneK` or `geneC` in `gene_counts.tsv`.

```bash
grep "^gene[KC]" gene_counts.tsv | wc -l
```

- Count how many sequences are in `random.fasta`.

```bash
grep -c "^>" random.fasta
```

---

### 4. **Text Processing with `sed`**
- Replace all instances of `gene` in `gene_counts.tsv` with `GENE`.

```bash
sed 's/gene/GENE/g' gene_counts.tsv | head
```

- Convert `data.tsv` into a CSV file.

```bash
sed 's/\t/,/g' data.tsv > data.csv
```

---

### 5. **Character-Level Transformations with `tr`**
- Convert all nucleotide sequences in `random.fasta` to lowercase.

```bash
tr 'A-Z' 'a-z' < random.fasta > random_lowercase.fasta
```

- Remove numbers from FASTA sequences.

```bash
tr -d '0-9' < random.fasta | head
```

---

### 6. **BONUS: Count Features in GTF File**
If you downloaded a GTF file from Ensembl:

```bash
wget ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

zcat Homo_sapiens.GRCh38.110.gtf.gz | grep -P '\tgene\t' | grep 'rRNA' | wc -l
zcat Homo_sapiens.GRCh38.110.gtf.gz | grep -P '\texon\t' | grep 'tRNA' | wc -l
```

---

### ✅ Submission (Optional)
Summarize your findings or output in a `report.txt` file using `echo` and redirection.

```bash
echo "Number of gene rows: $(wc -l < gene_counts.tsv)" > report.txt
echo "Number of sequences: $(grep -c '^>' random.fasta)" >> report.txt
```

---

**📂 Files Used:**
- `data.tsv`
- `gene_counts.tsv`
- `genes.txt`
- `lengths.txt`
- `random.fasta`

---

📌 Tip: Try modifying the scripts or extending them. For instance, filter genes with expression above 200, or find sequences longer than 10 bp using `awk`.


---

## 🧪 Advanced Practice Tasks

### Task 6: Filter Genes with High Expression in Multiple Samples

Using `awk`, filter rows in `gene_counts.tsv` where:
- **Sample3** > 250
- **Sample6** > 250

<details>
<summary>Click to view code</summary>

<pre><code class="language-bash">
awk '$4 > 250 && $7 > 250' destination/gene_counts.tsv
</code></pre>
</details>

---

### Task 7: Convert FASTA sequences to lowercase and count base composition

Use `tr` to convert the sequences in `random.fasta` to lowercase, and then use `grep` and `wc` to count how many 'a', 't', 'g', and 'c' characters are present.

<details>
<summary>Click to view code</summary>

<pre><code class="language-bash">
tr 'A-Z' 'a-z' < destination/random.fasta | grep -v "^>" | tee seq_lc.txt
grep -o "a" seq_lc.txt | wc -l
grep -o "t" seq_lc.txt | wc -l
grep -o "g" seq_lc.txt | wc -l
grep -o "c" seq_lc.txt | wc -l
</code></pre>
</details>

---

### Task 8: Merge and Sort Gene Data

Combine `data.tsv` and `genes.tsv` and sort the result based on the first column alphabetically.

<details>
<summary>Click to view code</summary>

<pre><code class="language-bash">
cat destination/data.tsv destination/genes.tsv > merged.tsv
sort merged.tsv
</code></pre>
</details>

---

### Task 9: Replace Underscores in Sequence Headers

In `random.fasta`, use `sed` to replace underscores `_` with spaces in the headers (lines starting with `>`).

<details>
<summary>Click to view code</summary>

<pre><code class="language-bash">
sed -E 's/^>([^_]+)_([^_]+)/>\1 \2/' destination/random.fasta | head
</code></pre>
</details>

---

### Task 10: Compress Multiple Files and Verify Sizes

Compress `data.tsv`, `genes.tsv`, and `gene_counts.tsv` using `gzip`, and compare the file sizes before and after.

<details>
<summary>Click to view code</summary>

<pre><code class="language-bash">
ls -lh destination/data.tsv destination/genes.tsv destination/gene_counts.tsv
gzip -k destination/*.tsv
ls -lh destination/*.gz
</code></pre>
</details>
