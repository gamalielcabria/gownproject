# 🧬 Bioinformatics Terminal Activity: Parsing Gene Data

This exercise helps students practice using `grep`, `cut`, and `awk` to parse and analyze tabular gene data using a simulated `genes.tsv` file.

---

## 🗂️ Data Setup

Run the following block in your terminal to create a test file with gene data:

```bash
cat > genes.tsv <<EOF
GeneID	Species	Function	Length_bp
EOF

for i in $(seq -w 1 100); do
  gene="G$(printf '%04d' $i)"
  species=$(shuf -n 1 -e Homo_sapiens Escherichia_coli Arabidopsis_thaliana Saccharomyces_cerevisiae Mus_musculus)
  function=$(shuf -n 1 -e kinase transporter hydrolase ligase polymerase synthetase)
  length=$((RANDOM % 1500 + 500))
  echo -e "${gene}	${species}	${function}	${length}" >> genes.tsv
done
```

---

## 🎯 Student Tasks

### 1. **Find all genes from *Homo sapiens***
```bash
grep "Homo_sapiens" genes.tsv
```

### 2. **Extract only the GeneID and Function columns**
```bash
cut -f1,3 genes.tsv
```

### 3. **List genes longer than 1000 bp**
```bash
awk 'NR==1 || $4 > 1000' genes.tsv
```

### 4. **Functions of *E. coli* genes longer than 850 bp**
```bash
grep "Escherichia_coli" genes.tsv | awk '$4 > 850 { print $3 }'
```

### 5. **Count how many genes are "kinase"**
```bash
awk '$3 == "kinase" { count++ } END { print count }' genes.tsv
```

---

## 💡 Bonus Tasks

- Print the list of **unique functions**:
  ```bash
  cut -f3 genes.tsv | sort | uniq
  ```

- Count how many genes each **species** has:
  ```bash
  awk 'NR>1 { count[$2]++ } END { for (s in count) print s, count[s] }' genes.tsv
  ```

- Find the **average gene length** per species:
  ```bash
  awk 'NR>1 { sum[$2]+=$4; count[$2]++ } END { for (s in sum) print s, sum[s]/count[s] }' genes.tsv
  ```

---

## ✅ Skills Tested

- Pattern matching with `grep`
- Column extraction with `cut`
- Conditional filtering and aggregation with `awk`
- Combining UNIX text processing tools
