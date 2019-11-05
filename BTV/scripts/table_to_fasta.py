#!/usr/bin/env python

# script to take "TABLE OF TREE SEGS" from John White
# and collate sequences from each seg into fasta files
# 10.07.2019

# sometimes the spreadsheet gives an NCBI# to download
# other times it gives an isolate name which I'll get from my current files
# will do one segment at a time

import sys, os
import argparse
from Bio import Entrez
from Bio import SeqIO
import csv

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("take TABLE OF TREE SEGS from John in csv format and collate"
                                 " sequences into fasta files")

parser.add_argument('-t', '--table', type = str,
                    help = "TABLE OF TREE SEGS from John. Must be in csv format.")
parser.add_argument('-s', '--seq_file', type = str,
                    help = "Fasta file containing already obtained sequences")
parser.add_argument('-m', '--meta_file', type = str,
                    help = "Meta data file for the obtained sequences")
parser.add_argument('-g', '--seg', type = str,
                    help = "Segment to analyse. Should be in the format 'SEG-4'")


args = parser.parse_args()

# copy function to get NCBI records from 'retrieve_genbank_records.py'

Entrez.email = "matthewjneave1@gmail.com"

def retrieve_ncbi_record(ncbi_id):
    print("retrieving {} from NCBI".format(ncbi_id))
    new_handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta")
    seq_record = SeqIO.read(new_handle, "fasta")
    return(seq_record)

# call the sequence file into an index for easy access

avail_seqs = SeqIO.index(args.seq_file, "fasta")

# make a meta data dictionary for the already obtained sequences

meta_dict = {}

with open(args.meta_file) as fl:
    header = next(fl).strip().split(",")
    for line in fl:
        line = line.strip()
        cols = line.split(",")
        strain = cols[0]
        meta_dict[strain] = cols[5:]

# read through spreadsheet line by line
# extracting appropriate sequence and metadata

records_to_write = []
meta_to_write = [["ID", "strain", "type", "year", "country", "accession"]]

with open(args.table) as t:
    reader = csv.reader(t)
    # only want header columns that refer to sequences
    # see next for loop
    header = next(reader)[4:]
    for line in reader:
        strain = line[0].strip()
        type = line[1].strip()
        year = line[2].strip()
        country = line[3].strip().replace("_", "/")
        for index, seq in enumerate(line[4:]):
            seq = seq.strip()
            # check the right segment is selected for this run
            # and that John wants this particular segment in the tree
            if seq != "" and header[index] == args.seg:
                # this is a sequence that we need to include
                ## 16.07.19 John has sent me a new spreadsheet that contains
                ## a new isolate name (replacing an old one)
                ## plus it has genbank numbers but these are not on NCBI
                ## need to hack my script to get this to work
                if strain == "DPP9244":
                    strain = "V9244"
                    new_name = "DPP9244"

                # first check if I have it already
                if strain in avail_seqs:
                    segment_number = int(args.seg.split("-")[1])
                    #print("{} detected in avail_seqs".format(strain))
                    records_to_write.append(avail_seqs[strain])
                    # if the record is V9244, need to include accession numbers from the spreadsheet
                    # these are not on NCBI yet
                    if strain == "V9244":
                        accession = seq
                    elif strain != seq:
                        accession = seq
                    else:
                        try:
                            accession = meta_dict[strain][segment_number]
                        except:
                            print("couldn't get accession number for {}, setting to unassigned".format(strain))
                            accession = "MZXXXXXX"

                    if strain == "V9244":
                        meta_to_write.append([strain, new_name, type, year, country, accession])
                    else:
                        meta_to_write.append([strain, strain, type, year, country, accession])
                else:
                    try:
                        record = retrieve_ncbi_record(seq)
                        records_to_write.append(record)
                        # note: genbank records get downloaded as fasta with the version number
                        # need to add this '.1' to the meta (e.g. KX578963.1)
                        meta_to_write.append([record.id, strain, type, year, country, record.id])
                    except Exception as e:
                        print(e)
                        print("{} unable to be collated".format(seq))

SeqIO.write(records_to_write, args.seg + ".fasta", "fasta")

with open(args.seg + ".meta", "w") as output:
    for meta in meta_to_write:
        output.write(",".join(meta) + "\n")












