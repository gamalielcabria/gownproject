#DADA2_exploratory.R
library(dada2)
library(tidyverse)

input <- "/home/glbcabria/Workbench/P0/qc" #Contains raw *fastq.gz 

# Taking in inputs for SLURM
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  n_cores <- 12
  cat("No number of cores provided. Defaulting to 12 cores.\n")
} else {
  n_cores <- as.numeric(args[1])
  if (is.na(n_cores) || n_cores < 1) {
    stop("Invalid number of cores.")
  }
}

# Listing inputs
forward_fq <- sort(list.files(input, pattern="_R1_paired.trim.fastq.gz", full.names = TRUE))
reverse_fq <- sort(list.files(input, pattern="_R2_paired.trim.fastq.gz", full.names = TRUE))
sample_names <- sapply(strsplit(basename(forward_fq), "_"), `[`, 1)

list_reads <- c(forward_fq, reverse_fq)

for (i in seq_along(forward_fq)){
    out_name <- paste0(input,'/',sample_names[i], ".qualplot.png")
    png(out_name, width=800, height = 1000)
    print(plotQualityProfile(c(forward_fq[i],reverse_fq[i])))
    dev.off()
}

# Filter and Trimming
filt_for_fq <- file.path(input, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
filt_rev_fq <- file.path(input, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))
names(filt_for_fq) <- sample_names
names(filt_rev_fq) <- sample_names

filtntrim <- filterAndTrim(
    forward_fq, filt_for_fq,
    reverse_fq, filt_rev_fq,
    truncLen = c(260,260),
    #maxReads = 100000, #Comment out on first run to see the distribution of reads
    maxN=0, maxEE=c(2,2), truncQ=2, 
    rm.phix=TRUE,
    compress=TRUE, 
    multithread=n_cores
)

# Error Rate
forward_err <- learnErrors(filt_for_fq, multithread = TRUE)
reverse_err <- learnErrors(filt_rev_fq, multithread = TRUE)

plot_errs_for <- plotErrors(forward_err, nominalQ = TRUE)
plot_errs_rev <- plotErrors(reverse_err, nominalQ = TRUE)
ggsave(plot = plot_errs_for, filename = paste0(input,'/plot_errs_for.png'))
ggsave(plot = plot_errs_rev, filename = paste0(input,'/plot_errs_rev.png'))

# Sample Inference
dada_for <- dada(filt_for_fq, forward_err, multithread = n_cores, pool = TRUE, DETECT_SINGLETONS=TRUE)
dada_rev <- dada(filt_rev_fq, reverse_err, multithread = n_cores, pool = TRUE, DETECT_SINGLETONS=TRUE)

dada_for[[1]]
denoisedF<-sapply(dada_for, function(x) sum(getUniques(x)))
denoisedR<-sapply(dada_rev, function(x) sum(getUniques(x)))

# Inspecting the dada-class object can show how many sequence variants from number of input sequences
merged_reads <- mergePairs(dada_for, filt_for_fq, dada_rev, filt_rev_fq, verbose=TRUE)

# Constructing sequence table
seqtab <- makeSequenceTable(merged_reads)
dim(seqtab)
merged_counts <- sapply(merged_reads, function(x) sum(x$abundance))
print(merged_counts)

# Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Tracking
getN <- function(x) sum(getUniques(x))
track <- cbind(filtntrim, sapply(dada_for, getN), sapply(dada_rev, getN), sapply(merged_reads, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names

final_track <- track %>% as.data.frame %>% 
mutate(Percent_Merged = merged / filtered ) %>% 
mutate(Percent_NoChimeras = nonchim / filtered )

write.csv(final_track, file = 'track.csv')
