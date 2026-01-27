#!/usr/bin/env python3
import sys

lens=[]
cur=0
for line in sys.stdin:
    if line.startswith(">"):
        if cur: lens.append(cur)
        cur=0
    else:
        cur += len(line.strip())
if cur: lens.append(cur)

lens.sort(reverse=True)
total=sum(lens)
half=total/2
run=0
for L in lens:
    run += L
    if run >= half:
        print(L)
        break
