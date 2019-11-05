#!/usr/bin/env python

# split genbank file into structural / non-structural for nextstrain
# 13.12.2018

import sys
from Bio import SeqIO

gb_file = sys.argv[1]

record = SeqIO.read(open(gb_file), "genbank")

non_struct = record[:5277]
struct = record[5277:]

SeqIO.write(non_struct, "KT280060_NS.gb", "genbank")
SeqIO.write(struct, "KT280060_S.gb", "genbank")
