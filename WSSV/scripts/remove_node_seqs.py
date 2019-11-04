#!/usr/bin/env python
# python 3

# in the gene translations, the ancestral 'node' sequences are included
# want to remove these and just give the real sequences to nextstrain

import sys
from Bio import SeqIO

trans_fasta = sys.argv[1]
output = open(sys.argv[2], "w")

for record in SeqIO.parse(trans_fasta, "fasta"):
    if not record.id.startswith("NODE"):
        SeqIO.write(record, output, "fasta")
