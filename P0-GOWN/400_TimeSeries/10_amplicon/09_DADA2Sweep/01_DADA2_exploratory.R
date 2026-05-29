####!/usr/bin/env Rscript
# DADA2_exploratory.R
# Usage:
#   Rscript DADA2_exploratory.R <INPUT_DIR> <OUTPUT_DIR> [TAX_TRAIN_FASTA] [SPECIES_FASTA] [TRUNC_F] [TRUNC_R]

suppressPackageStartupMessages({
  library(dada2)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript DADA2_exploratory.R <INPUT_DIR> <OUTPUT_DIR> [TAX_TRAIN_FASTA] [SPECIES_FASTA] [TRUNC_F] [TRUNC_R]")
}

input_dir  <- args[[1]]
output_dir <- args[[2]]

# Optional taxonomy files kept here only so script signature stays compatible
tax_dir <- "/work/ebg_lab/eb/AIWIP_GOWN_amplicons/001_scripts/DADA2/assignedTaxa"
default_tax_train_fa <- file.path(tax_dir, "silva_nr99_v138.2_toGenus_trainset.fa.gz")
default_species_fa   <- file.path(tax_dir, "silva_v138.2_assignSpecies.fa.gz")

tax_train_fa <- if (length(args) >= 3 && nzchar(args[[3]])) args[[3]] else default_tax_train_fa
species_fa   <- if (length(args) >= 4 && nzchar(args[[4]])) args[[4]] else default_species_fa

# New: optional truncLen values
truncF <- if (length(args) >= 5 && nzchar(args[[5]])) as.integer(args[[5]]) else 250L
truncR <- if (length(args) >= 6 && nzchar(args[[6]])) as.integer(args[[6]]) else 230L

if (is.na(truncF) || is.na(truncR) || truncF < 1 || truncR < 1) {
  stop("Invalid truncLen values: truncF=", truncF, " truncR=", truncR)
}

threads <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA_character_)
threads <- if (is.na(threads) || threads == "") 12L else as.integer(threads)
if (is.na(threads) || threads < 1) threads <- 12L

# Output folders
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(output_dir, "plots")
filt_dir  <- file.path(output_dir, "filtered")
rds_dir   <- file.path(output_dir, "rds")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(filt_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir,   recursive = TRUE, showWarnings = FALSE)

cat("========== DADA2 RUN ==========\n")
cat("INPUT_DIR       :", input_dir, "\n")
cat("OUTPUT_DIR      :", output_dir, "\n")
cat("THREADS         :", threads, "\n")
cat("TRUNC_F         :", truncF, "\n")
cat("TRUNC_R         :", truncR, "\n")
cat("RDS_DIR         :", rds_dir, "\n\n")

writeLines(
  c(
    paste("INPUT_DIR =", input_dir),
    paste("OUTPUT_DIR =", output_dir),
    paste("THREADS =", threads),
    paste("TRUNC_F =", truncF),
    paste("TRUNC_R =", truncR)
  ),
  con = file.path(output_dir, "run_parameters.txt")
)

# ---------------------------
# Find reads
# ---------------------------
fwd_pat <- "_R1_001\\.fastq\\.gz$"
rev_pat <- "_R2_001\\.fastq\\.gz$"

forward_fq <- sort(list.files(input_dir, pattern = fwd_pat, full.names = TRUE))
reverse_fq <- sort(list.files(input_dir, pattern = rev_pat, full.names = TRUE))

if (length(forward_fq) == 0 || length(reverse_fq) == 0) {
  stop("No FASTQ files detected in input_dir with expected patterns.")
}

sample_names <- sub(fwd_pat, "", basename(forward_fq))

# ---------------------------
# Filtering
# ---------------------------
cat("Filtering + trimming...\n")

filt_for_fq <- file.path(filt_dir, paste0(sample_names, "_F_filt.fastq.gz"))
filt_rev_fq <- file.path(filt_dir, paste0(sample_names, "_R_filt.fastq.gz"))
names(filt_for_fq) <- sample_names
names(filt_rev_fq) <- sample_names

filtntrim <- filterAndTrim(
  forward_fq, filt_for_fq,
  reverse_fq, filt_rev_fq,
  truncLen = c(truncF, truncR),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = threads
)

write.csv(filtntrim, file = file.path(output_dir, "01_filtering_summary.csv"))
saveRDS(filtntrim, file = file.path(rds_dir, "01_filtntrim.rds"))

# ---------------------------
# Error model
# ---------------------------
cat("Learning error models...\n")

forward_err <- learnErrors(filt_for_fq, multithread = threads)
reverse_err <- learnErrors(filt_rev_fq, multithread = threads)

saveRDS(forward_err, file = file.path(rds_dir, "02_forward_err.rds"))
saveRDS(reverse_err, file = file.path(rds_dir, "02_reverse_err.rds"))

pdf(file.path(plots_dir, "02_error_models.pdf"))
plotErrors(forward_err, nominalQ = TRUE)
plotErrors(reverse_err, nominalQ = TRUE)
dev.off()

# ---------------------------
# DADA inference
# ---------------------------
cat("Running DADA inference...\n")

derepF <- derepFastq(filt_for_fq)
derepR <- derepFastq(filt_rev_fq)
names(derepF) <- sample_names
names(derepR) <- sample_names

dada_for <- dada(filt_for_fq, err = forward_err, multithread = threads, pool = "pseudo")
dada_rev <- dada(filt_rev_fq, err = reverse_err, multithread = threads, pool = "pseudo")

saveRDS(dada_for, file = file.path(rds_dir, "04_dada_for.rds"))
saveRDS(dada_rev, file = file.path(rds_dir, "04_dada_rev.rds"))

# ---------------------------
# Merge pairs
# ---------------------------
cat("Merging pairs...\n")

merged_reads <- mergePairs(dada_for, derepF, dada_rev, derepR, verbose = TRUE)
saveRDS(merged_reads, file = file.path(rds_dir, "05_merged_reads.rds"))

merged_counts <- sapply(merged_reads, function(x) sum(x$abundance))
write.csv(
  data.frame(Sample = names(merged_counts), merged = as.integer(merged_counts)),
  file = file.path(output_dir, "05_merged_counts.csv"),
  row.names = FALSE
)

# ---------------------------
# Sequence table
# ---------------------------
cat("Building sequence table...\n")

seqtab <- makeSequenceTable(merged_reads)
saveRDS(seqtab, file = file.path(rds_dir, "06_seqtab_raw.rds"))

# ---------------------------
# Chimera removal
# ---------------------------
cat("Removing chimeras...\n")

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = threads,
  verbose = TRUE
)
saveRDS(seqtab.nochim, file = file.path(rds_dir, "07_seqtab_nochim.rds"))

# ---------------------------
# Tracking
# ---------------------------
cat("Creating tracking table...\n")

getN <- function(x) sum(getUniques(x))

track <- data.frame(
  Sample    = sample_names,
  input     = filtntrim[, 1],
  filtered  = filtntrim[, 2],
  denoisedF = sapply(dada_for, getN),
  denoisedR = sapply(dada_rev, getN),
  merged    = sapply(merged_reads, getN),
  nonchim   = rowSums(seqtab.nochim),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Percent_Merged     = if_else(filtered > 0, merged / filtered, NA_real_),
    Percent_NoChimeras = if_else(filtered > 0, nonchim / filtered, NA_real_)
  )

write.csv(track, file = file.path(output_dir, "06_track.csv"), row.names = FALSE)
saveRDS(track, file = file.path(rds_dir, "08_track.rds"))

# ---------------------------
# OTU/ASV table
# ---------------------------
cat("Writing OTU/ASV table...\n")

otu_table <- as.data.frame(seqtab.nochim)
write.csv(otu_table, file = file.path(output_dir, "08_otu_table.csv"), row.names = TRUE)
saveRDS(otu_table, file = file.path(rds_dir, "10_otu_table.rds"))

cat("DADA2 pipeline completed.\n")
cat("Outputs written to: ", output_dir, "\n")
cat("RDS written to    : ", rds_dir, "\n")
