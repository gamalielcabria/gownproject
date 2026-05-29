#!/bin/bash

echo -e "combination\tmean_percent_merged\tmean_percent_nochimeras\tmedian_percent_merged\tmedian_percent_nochimeras"

for d in truncF*_truncR*; do
    f="$d/06_track.csv"
    [[ -f "$f" ]] || continue

    awk -F',' -v combo="$d" '
    NR==1 {next}

    {
        gsub(/"/, "", $8)
        gsub(/"/, "", $9)

        if ($8 != "" && $8 != "NA") {
            n8++
            sum8 += $8
            a8[n8] = $8
        }

        if ($9 != "" && $9 != "NA") {
            n9++
            sum9 += $9
            a9[n9] = $9
        }
    }

    END {
        if (n8 > 0 && n9 > 0) {

            for (i=1; i<=n8; i++) {
                for (j=i+1; j<=n8; j++) {
                    if (a8[i] > a8[j]) {
                        tmp = a8[i]
                        a8[i] = a8[j]
                        a8[j] = tmp
                    }
                }
            }

            for (i=1; i<=n9; i++) {
                for (j=i+1; j<=n9; j++) {
                    if (a9[i] > a9[j]) {
                        tmp = a9[i]
                        a9[i] = a9[j]
                        a9[j] = tmp
                    }
                }
            }

            mean8 = sum8 / n8
            mean9 = sum9 / n9

            if (n8 % 2) med8 = a8[(n8+1)/2]
            else        med8 = (a8[n8/2] + a8[n8/2+1]) / 2

            if (n9 % 2) med9 = a9[(n9+1)/2]
            else        med9 = (a9[n9/2] + a9[n9/2+1]) / 2

            printf "%s\t%.6f\t%.6f\t%.6f\t%.6f\n", combo, mean8, mean9, med8, med9
        }
    }
    ' "$f"
done
