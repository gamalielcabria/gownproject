#!/usr/bin/env bash
set -euo pipefail

CSV="/work/ebg_lab/Instrument_data/Sequences/MGenomes/run251210/GSC-2651_sequence-name-index.csv"   # <- change this
OUTDIR="/work/ebg_lab/eb/AIWIP_GOWN_metagenome/100_denitrificationindex/00_rawdata"           # <- change if you want
INPUT_DIR="/work/ebg_lab/Instrument_data/Sequences/MGenomes/run251210/"  # <- change this
mkdir -p "$OUTDIR"

# Read CSV (skip header). Strip Windows CRs if present.
tail -n +2 "$CSV" | tr -d '\r' | \
while IFS=',' read -r ShortCode Name NEB_Code ForID ForSeq RevID RevSeq; do
  # Skip empty/short lines
  [[ -z "${Name:-}" || -z "${ForSeq:-}" || -z "${RevSeq:-}" ]] && continue

  # Drop trailing Ns from ForSeq (e.g., CAAGGTACNNNN... -> CAAGGTAC)
  ForNoN="${ForSeq%%N*}"
  barcode="${ForNoN}-${RevSeq}"

  r1_glob="${INPUT_DIR}/*_*_1_${barcode}_150bp.concat.fastq.gz"
  r2_glob="${INPUT_DIR}/*_*_2_${barcode}_150bp.concat.fastq.gz"

  # Collect matches using compgen (won't leave literal globs when no match)
  mapfile -t r1 < <(compgen -G "$r1_glob" || true)
  mapfile -t r2 < <(compgen -G "$r2_glob" || true)

  if [[ ${#r1[@]} -ne 1 || ${#r2[@]} -ne 1 ]]; then
    echo "WARN: ${Name} (${barcode}) -> R1 matches=${#r1[@]} R2 matches=${#r2[@]} (skipping)" >&2
    echo "      R1 pattern: $r1_glob" >&2
    echo "      R2 pattern: $r2_glob" >&2
    continue
  fi

  # FIX 2: explicitly require that the FASTQs exist before linking
  if [[ ! -e "${r1[0]}" || ! -e "${r2[0]}" ]]; then
    echo "WARN: Missing FASTQ(s) for ${Name} (${barcode}) (skipping)" >&2
    echo "      R1: ${r1[0]}" >&2
    echo "      R2: ${r2[0]}" >&2
    continue
  fi

  link_r1="${OUTDIR}/${Name}__R1_${barcode}.fastq.gz"
  link_r2="${OUTDIR}/${Name}__R2_${barcode}.fastq.gz"

  # Create/overwrite symlinks
  ln -sfn "$(readlink -f "${r1[0]}")" "$link_r1"
  ln -sfn "$(readlink -f "${r2[0]}")" "$link_r2"

  echo "OK: ${Name} (${barcode})"
done