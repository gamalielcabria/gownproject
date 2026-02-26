#!/usr/bin/env bash
set -euo pipefail

SRCDIR="."
OUTFILE="fastq_readcount_table.tsv"
RESUMEFILE=""
OVERWRITE=0

usage() {
  cat <<EOF
Usage:
  $0 [-s DIR] [-o OUT.tsv] [-r RESUME.tsv] [--overwrite]

Options:
  -s, --src        Source directory of fastq.gz files (default: .)
  -o, --out        Output TSV file (default: fastq_readcount_table.tsv)
  -r, --resume     Optional existing TSV to resume from (skip already-counted R1s)
      --overwrite  Overwrite OUTFILE if it already exists
EOF
}

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--src) SRCDIR="$2"; shift 2 ;;
    -o|--out) OUTFILE="$2"; shift 2 ;;
    -r|--resume) RESUMEFILE="$2"; shift 2 ;;
    --overwrite) OVERWRITE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

header=$'R1_file\tR2_file\tR1_reads\tR2_reads\tMatch'

if [[ -e "$OUTFILE" && "$OVERWRITE" -ne 1 ]]; then
  echo "ERROR: $OUTFILE exists. Use --overwrite or choose a different -o output." >&2
  exit 1
fi

# Start output with header
printf "%s\n" "$header" > "$OUTFILE"

declare -A done

# If a resume file is provided and non-empty, copy rows + load completed R1s
if [[ -n "$RESUMEFILE" ]]; then
  if [[ ! -s "$RESUMEFILE" ]]; then
    echo "ERROR: Resume file provided but not found or empty: $RESUMEFILE" >&2
    exit 1
  fi

  # Append existing rows (skip header line)
  tail -n +2 "$RESUMEFILE" >> "$OUTFILE"

  # Load completed R1s into hash
  while IFS=$'\t' read -r r1name _; do
    [[ -n "$r1name" ]] && done["$r1name"]=1
  done < <(awk -F'\t' 'NR>1 && $1!="" {print $1 "\t" $2}' "$RESUMEFILE")
fi

shopt -s nullglob

for r1 in "$SRCDIR"/*_R1_*.fastq.gz; do
  r1base="$(basename "$r1")"

  # Skip if already counted (only matters when resuming)
  if [[ -n "${done[$r1base]+x}" ]]; then
    continue
  fi

  r2="${r1/_R1_/_R2_}"
  r2base="$(basename "$r2")"

  if [[ ! -f "$r2" ]]; then
    echo "WARN: Missing R2 for $r1base (skipping)" >&2
    continue
  fi

  echo "Processing: $r1base"

  r1_reads=$(gzip -dc "$r1" | awk 'END{print NR/4}')
  r2_reads=$(gzip -dc "$r2" | awk 'END{print NR/4}')

  if [[ "$r1_reads" -eq "$r2_reads" ]]; then
    status="OK"
  else
    status="MISMATCH"
    echo "WARNING: Read count mismatch for $r1base" >&2
  fi

  printf "%s\t%s\t%s\t%s\t%s\n" "$r1base" "$r2base" "$r1_reads" "$r2_reads" "$status" >> "$OUTFILE"
done

echo "Done. Output written to $OUTFILE"