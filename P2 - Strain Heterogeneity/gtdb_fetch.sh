#!/bin/bash

workdir=/work/ebg_lab/referenceDatabases/GTDB_R214
destination=/work/ebg_lab/eb/strain_heterogeneity/GTDB
version="`cat metadata.txt | grep '^VERSION_DATA' | cut -d'=' -f2`"

# Generate list of all genomes in the database

find "${workdir}" -name "*.fna.gz" | while read filepath; do
    echo -e "${filepath}\t$(basename "${filepath}" .fna.gz)"
done > genome_paths.list

# Make destination directory
mkdir -p destination

# Transfer the files 
while IFS=$'\t' read -r fullpath basename; do
    cp "$fullpath" "$destination/${basename}.fna.gz"
done < genome_paths.list
