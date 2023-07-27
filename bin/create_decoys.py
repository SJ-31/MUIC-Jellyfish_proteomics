#!/usr/bin/env python
import sys
from Bio import SeqIO

decoy: str = sys.argv[2]
for record in SeqIO.parse(sys.argv[1], "fasta"):
    with open(decoy, "a") as a:
        a.write(">decoy_" + record.id[1:] + "\n")
        a.write(str(record.seq[::-1]) + "\n")
