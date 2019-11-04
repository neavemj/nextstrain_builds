#!/usr/bin/env python

# split a multi-genbank file into fasta, plus associated metadata
# Matthew J. Neave 05.12.2018

import sys
from Bio import SeqIO # requires biopython

gb_file = sys.argv[1]
meta_file = open(sys.argv[2], "w")

for record in SeqIO.parse(open(gb_file), "genbank"):
    try:
        country = record.features[0].qualifiers["country"][0]
    except:
        print("** country not found for: {}. Assigning to unknown".format(record.name))
        country = "Unknown"
    try:
        collection_date = record.features[0].qualifiers["collection_date"][0]
    except:
        print("** collection date not found for: {}. Assigning to unknown.".format(record.name))
        collection_date = "Unknown"
    try:
        host = record.features[0].qualifiers["host"][0]
    except:
        print("** host not found for: {}. Assigning to unknown".format(record.name))
        host = "Unknown"

    authors = ",".join(record.annotations['references'][0].authors.split(",")[0:6]) + " et al"
    title = record.annotations['references'][0].title
    journal = record.annotations['references'][0].journal

    # my other meta file has these headers
    # strain	location	country	state	source	date	specimen	host	platform authors title journal

    meta_file.write("\t".join([record.name, "Unknown", country, "Unknown", "NCBI", collection_date, "Unknown", host,
                               "Unknown", authors, title, journal]) + "\n")

    fasta_file = open(record.name + ".fasta", "w")
    fasta_file.write(">" + record.name + "\n" + str(record.seq) + "\n")