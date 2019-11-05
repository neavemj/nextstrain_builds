#!/usr/bin/env python

# Add newly downloaded genbank info from exotic segs to current BTV build
# also update the lat long file if required

import sys
import argparse
from Bio import SeqIO
from datetime import datetime

parser = argparse.ArgumentParser("update a Nextstrain build with new exotic NCBI sequences\n"
                                 "meta and fasta files will be output for each segment")

parser.add_argument('-g', '--genbank_files', type = str,
        nargs = "+", help = "genbank files to convert to nextstrain")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check required arguments are provided

if args.genbank_files is None:
    print("\n** required input missing\n"
          "** genbank files are required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# function to parse downloaded genbank records

def parse_genbank(fl_name):
    seq_record = SeqIO.read(fl_name, "genbank")
    record_dict = {}

    # extracting several bits of information from different parts of the record
    # some are absolutely required for nextstrain others are optional
    required_qualifiers = ["strain", "segment", "collection_date"]
    optional_qualifiers = ["host", "serotype", "country", "lat_lon"]

    for req in required_qualifiers:
        try:
            info = seq_record.features[0].qualifiers[req][0]
            record_dict[req] = info
        except KeyError:
            print("{} not found in record {}.".format(req, fl_name))
            print("This is a required field - record aborted. \n")
            return(None)

    for opt in optional_qualifiers:
        try:
            info = seq_record.features[0].qualifiers[opt][0]
            record_dict[opt] = info
        except KeyError:
            print("{} not found in record {}.".format(req, fl_name))
            print("This is optional - record will be written with 'Unknown' in this field. \n")
            record_dict[opt] = "Unknown"

    record_dict["authors"] = ",".join(seq_record.annotations['references'][0].authors.split(",")[0:2]) + " et al"
    record_dict["title"] = seq_record.annotations['references'][0].title
    record_dict["journal"] = seq_record.annotations['references'][0].journal

    record_dict["sequence"] = str(seq_record.seq)

    return(record_dict)

# some bits of info in the genbank require some additional processing

def additional_processing(raw_record_dict):
    if raw_record_dict["country"] != "Unknown":
        country = raw_record_dict["country"].split(":")[0].lower()
        state = raw_record_dict["country"].split(":")[1].split(",")[1].strip()
        town = raw_record_dict["country"].split(":")[1].split(",")[0].strip().lower().replace(" ", "_")
    else:
        state, town = "Unknown"

    raw_record_dict["country"] = country
    raw_record_dict["state"] = state
    raw_record_dict["town"] = town
    raw_record_dict["place"] = town

    # also want the serotype to be a str, not int
    if raw_record_dict["serotype"] != "Unknown":
        raw_record_dict["serotype"] = "type_" + raw_record_dict["serotype"]

    # also want to convert the date properly
    # useful info available here:
    # https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior
    date = raw_record_dict["collection_date"]
    #python_date = datetime.strptime(date, "%b-%Y")
    python_date = datetime.strptime(date, "%d-%b-%Y")
    formatted_date = str(python_date.year) + "-" + str(python_date.month) + "-XX"
    raw_record_dict["collection_date"] = formatted_date

    return(raw_record_dict)

# extract all required info from the genbank files

genbank_dict = {}

for genbank_file in args.genbank_files:
    original_dict = parse_genbank(genbank_file)
    if original_dict is not None:
        genbank_dict[genbank_file] = additional_processing(original_dict)

# write out the lat longs with timestamp

today = datetime.today().strftime("%d-%m-%Y")

# now for the fasta and meta files
# using the 'a' parameter in the open function to append to a file rather than overwrite

for key in genbank_dict:
    meta_fl = open("btv_segment" + genbank_dict[key]["segment"] + "_" + today + "_meta.tsv", "a")
    fasta_fl = open("btv_segment" + genbank_dict[key]["segment"] + "_" + today + ".fasta", "a")
    cur = genbank_dict[key]

    meta_fl.write("\t".join([cur["strain"], "BTV", key.split(".")[0], cur["segment"], cur["serotype"], cur["state"], cur["collection_date"], cur["host"], cur["lat_lon"], cur["authors"], cur["title"], cur["journal"], cur["place"], cur["country"]]) + "\n")

    fasta_fl.write(">" + cur["strain"] + "\n" + cur["sequence"] + "\n")
