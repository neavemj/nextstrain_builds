#!/usr/bin/env python

# Add newly downloaded genbank info to current BTV build
# also update the lat long file if required

import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser("update a Nextstrain build with new NCBI sequences\n")

parser.add_argument('-g', '--genbank_files', type = str,
        nargs = "+", help = "genbank files to add (space separated)")
parser.add_argument('-l', '--current_lat_longs', type = str,
        nargs = "?", help = "the lat_long file for the existing build (new lat longs will be added if required)")
parser.add_argument('-f', '--current_fastas', type = str,
        nargs = "?", help = "the fasta files for the existing build (new sequences will be added)")
parser.add_argument('-m', '--current_metas', type = str,
        nargs = "?", help = "the meta files for the existing build (new meta data will be added)")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check required arguments are provided

if args.genbank_files is None or \
   args.current_lat_longs is None:
    print("\n** required input missing\n"
          "** a file containing wanted accession numbers is required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

acc_handle = args.acc_file

# function to parse downloaded genbank records

def parse_genbank(fl_name):
    seq_record = SeqIO.read(fl_name, "genbank")
    record_dict = {}

    # extracting several bits of information from different parts of the record
    # some are absolutely required for nextstrain others are optional
    required_qualifiers = ["isolate", "segment", "collection_date", "lat_lon"]
    optional_qualifiers = ["host", "serotype", "country"]

    for req in required_qualifiers:
        try:
            info = seq_record.features[0].qualifiers[req][0]
            record_dict[req] = info
        except KeyError:
            print("{} not found in record {}.".format(req, fl_name))
            print("This is a required field - record aborted. \n")

    for opt in optional_qualifiers:
        try:
            info = seq_record.features[0].qualifiers[opt][0]
            record_dict[opt] = info
        except KeyError:
            print("{} not found in record {}.".format(req, fl_name))
            print("This is an optional - record will be written with 'unknown' in this field. \n")
            record_dict[opt] = "Unknown"

    record_dict["authors"] = ",".join(seq_record.annotations['references'][0].authors.split(",")[0:3]) + " et al"
    record_dict["title"] = seq_record.annotations['references'][0].title
    record_dict["journal"] = seq_record.annotations['references'][0].journal

    return(record_dict)

feat = parse_genbank(sys.argv[1])
print(feat)
