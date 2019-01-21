#!/usr/bin/env python
# python 3

# for the subsetted WSSV builds, I also need to subset the gb file
# otherwise the alignment is weird and heaps of Ns get inserted

import sys
from Bio import SeqIO

orig_gb = sys.argv[1]
output = open(sys.argv[3], "w")

for record in SeqIO.parse(orig_gb, "genbank"):
    for feat in record.features:
        try:
            gene = feat.qualifiers["gene"]
            if gene[0] == sys.argv[2]:
                location = feat.location
                start = feat.location.start.position
                end = feat.location.end.position
                sub_record = record[start:end]
                SeqIO.write(sub_record, output, "genbank")
        except:
            pass
