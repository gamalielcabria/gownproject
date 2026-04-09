#!/bin/bash

for i in *taxa.txt; do
  taxa="$i"
  taxa_relabund="${taxa%.txt}.relabund.txt"

  awk '
  BEGIN{FS=OFS="\t"}
  NR==1 {print $0, "GW", "relabund_pct"; next}
  {
    samp = $1
    cov  = $2
    tax  = $3

    # extract GW ID from sample
    if (match(samp, /GW[0-9]+/)) {
      gw = substr(samp, RSTART, RLENGTH)
    } else {
      gw = "NA"
    }

    row[NR]    = $0
    sample[NR] = samp
    covar[NR]  = cov
    gw_id[NR]  = gw

    if (tax == "Root; d__Archaea" || tax == "Root; d__Bacteria") {
      total[samp] += cov
    }
  }
  END {
    for (i=2; i<=NR; i++) {
      if (total[sample[i]] > 0)
        print row[i], gw_id[i], 100 * covar[i] / total[sample[i]]
      else
        print row[i], gw_id[i], 0
    }
  }
  ' "$taxa" > "$taxa_relabund"

done
