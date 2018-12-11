#!/usr/bin/env python

# Robyn's fasta and metadata will need a bit of tidying before Nextstrain
# Matthew J. Neave 11.12.2018

# this script does the following:
# splits the file into RHDV genotype - RHDV1, RHDV2 and RCV
# changes U to T
# amends strain name (AUS/TAS/LEG-2) rather than NCBI accession as 'strain' ID

import sys
from Bio import SeqIO

fasta_file = sys.argv[1]
meta_file = sys.argv[2]

variant_convert = {"GI.1a": "RHDV1", "GI.1c": "RHDV1", "GI.4eP_GI.1a": "RHDV1", "GI.1bP_GI.2": "RHDV2", "GI.4eP_GI.2":
    "RCV"}

# need to write meta data separately for each variant

RHDV1_fasta = open("RHDV_RHDV1.fasta", "w")
RHDV1_meta = open("RHDV_RHDV1.meta.tsv", "w")
RHDV2_fasta = open("RHDV_RHDV2.fasta", "w")
RHDV2_meta = open("RHDV_RHDV2.meta.tsv", "w")
RCV_fasta = open("RHDV_RCV.fasta", "w")
RCV_meta = open("RHDV_RCV.meta.tsv", "w")

for f in [RHDV1_meta, RHDV2_meta, RCV_meta]:
    f.write("\t".join(["strain", "strain_short", "accession", "variant_long", "variant", "date", "country", "state", "authors",
                       "title", "host"]) + "\n")

# want a accession to strain dict for the fasta writing

acc_to_strain = {}
strain_to_variant = {}
strain_list = []

with open(meta_file) as f:
    next(f)
    for line in f:
        line = line.strip()
        cols = line.split("\t")
        strain = cols[0].strip()
        strain_list.append(strain)
        strain_short = cols[1]
        acc = cols[3].strip()
        variant_long = cols[4].strip()
        variant = variant_convert[variant_long]
        strain_to_variant[strain] = variant
        date = cols[5]
        country = cols[6]
        state = cols[7]
        authors = cols[8]
        title = cols[9]
        host = cols[10]
        lat = cols[11]
        long = cols[12]
        lat_long = lat + " " + long

        acc_to_strain[acc] = strain
        meta_list = [strain, strain_short, acc, variant_long, variant, date, country, state, authors, title, host,
                     lat_long]
        meta_list = [item.strip() for item in meta_list]
        if variant == "RHDV1":
            RHDV1_meta.write("\t".join(meta_list) + "\n")
        elif variant == "RHDV2":
            RHDV2_meta.write("\t".join(meta_list) + "\n")
        elif variant == "RCV":
            RCV_meta.write("\t".join(meta_list) + "\n")
        else:
            print("warning unknown variant in meta file: {}".format(variant))

# now go through the fasta file and split into variants with correct strain names

count = 0

def check_header(header):

    # sometimes the header in the fasta file contains an additional 2 digits at the end
    header_sub = "/".join(header.split("/")[:-1])

    # a couple of the strains have the accession number / strain
    possible_acc = header.split("/")[0]

    # some of the 'short' strain names have case differences
    header_upper = header_sub.upper()

    if header_sub in strain_list:
        return(header_sub)

    # check if the header is an NCBI accession
    elif header in acc_to_strain:
        return(acc_to_strain[header])

    # a couple of the strains have the accession number / strain
    elif possible_acc in acc_to_strain:
        return(acc_to_strain[possible_acc])

    elif header_upper in strain_list:
        return(header_upper)

    else:
        print("could not find header: {}, {}".format(header, header_sub))
        return("not found")


with open(fasta_file) as f:
    for record in SeqIO.parse(f, "fasta"):
        header = record.id.lstrip("O.cun*/").lstrip("L.eur*/").split(" ")[0]
        proper_header = check_header(header)
        record.id = proper_header
        record.description = proper_header
        # apparently biopython has a module for turning U back to T
        record.seq = record.seq.back_transcribe()
        variant = strain_to_variant[proper_header]
        if variant == "RHDV1":
            SeqIO.write(record, RHDV1_fasta, "fasta")
        elif variant == "RHDV2":
            SeqIO.write(record, RHDV2_fasta, "fasta")
        elif variant == "RCV":
            SeqIO.write(record, RCV_fasta, "fasta")
        else:
            print("warning unknown variant in fasta file: {}".format(variant))


