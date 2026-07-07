#!/usr/bin/env bash
set -euo pipefail

usage() {
cat << EOF
Edit FASTA headers by appending the filename.

USAGE
    $(basename "$0") -i INPUT_DIR [-o OUTPUT_DIR]

OPTIONS
    -i   Directory containing FASTA files (*.fasta)      [required]
    -o   Output directory for renamed FASTA files        [default: INPUT_DIR]
    -h   Show this help message

DESCRIPTION
    For each FASTA file, the script modifies headers so that:
        >GCA_013839515_1        398 bp          ami

    becomes
        >GCA_013839515_1--FILENAME

    where FILENAME is the FASTA file name without extension.
EXAMPLE
    $(basename "$0") -i hmms_nirs
    $(basename "$0") -i hmms_nirs -o renamed_fastas
EOF
}

INPUT_DIR=""
OUTPUT_DIR=""

while getopts "i:o:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [[ -z "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory required."
    usage
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR="$INPUT_DIR"
fi

mkdir -p "$OUTPUT_DIR"

shopt -s nullglob
files=("$INPUT_DIR"/*.fasta)

if [[ ${#files[@]} -eq 0 ]]; then
    echo "No FASTA files found in $INPUT_DIR"
    exit 1
fi

for fasta in "${files[@]}"; do
    base=$(basename "$fasta")
    name="${base%.fasta}"

    awk -v n="$name" '
    /^>/ {
        split($0, a, /[ \t]+/)
        print ">" a[1] "--" n
        next
    }
    {print}
    ' "$fasta" > "$OUTPUT_DIR/${name}_renamed.fasta"

    echo "Processed $base -> ${name}_renamed.fasta"
done

echo "Done."
