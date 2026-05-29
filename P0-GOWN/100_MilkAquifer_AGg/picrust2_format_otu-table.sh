#!/usr/bin/bash

IN="$1"

awk -F',' '
BEGIN{
    OFS="\t"
}

NR==1{
    # store ASV sequences
    for(i=2;i<=NF;i++){
        gsub(/^"|"$/, "", $i)
        seq[i]="ASV" i-1
        fasta[i]=$i
    }
    next
}

{
    gsub(/^"|"$/, "", $1)
    samples[NR]=$1

    for(i=2;i<=NF;i++){
        counts[i,NR]=$i
    }
}

END{
    # TSV header
    printf "ASV"
    for(r=2;r<=NR;r++){
        printf OFS samples[r]
    }
    printf "\n"

    # abundance table
    for(i=2;i<=length(seq)+1;i++){
        printf seq[i]
        for(r=2;r<=NR;r++){
            printf OFS counts[i,r]
        }
        printf "\n"
    }

    # FASTA output
    for(i=2;i<=length(seq)+1;i++){
        print ">" seq[i] "\n" fasta[i] > "rep_seqs.fasta"
    }
}
' $IN
