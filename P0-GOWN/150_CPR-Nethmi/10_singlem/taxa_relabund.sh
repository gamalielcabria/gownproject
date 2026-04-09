#!/bin/bash

for i in `ls *taxa.txt`;
do
taxa="$i"
taxa_relabund="${taxa%.txt}.relabund.txt"

awk '
BEGIN{FS=OFS="\t"}
NR==1 {print $0, "relabund_pct"; next}
{
  samp = $1
  cov  = $2
  tax  = $3

  row[NR]    = $0
  sample[NR] = samp
  covar[NR]  = cov

  if (tax == "Root; d__Archaea" || tax == "Root; d__Bacteria") {
    total[samp] += cov
  }
}
END {
  for (i=2; i<=NR; i++) {
    if (total[sample[i]] > 0) print row[i], 100 * covar[i] / total[sample[i]]
    else                      print row[i], 0
  }
}
' "$taxa" > "$taxa_relabund"

done
