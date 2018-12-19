#!/usr/bin/env python

# Add newly downloaded genbank info to current BTV build
# also update the lat long file if required

import sys
import argparse
from Bio import SeqIO
from datetime import datetime

parser = argparse.ArgumentParser("update a Nextstrain build with new NCBI sequences\n")

parser.add_argument('-g', '--genbank_files', type = str,
        nargs = "+", help = "genbank files to add (space separated)")
parser.add_argument('-l', '--current_lat_longs', type = str,
        nargs = "?", help = "the lat_long file for the existing build (new lat longs will be added if required)")
parser.add_argument('-f', '--current_fastas', type = str,
        nargs = "+", help = "the fasta files for the existing build (space separated)")
parser.add_argument('-m', '--current_metas', type = str,
        nargs = "+", help = "the meta files for the existing build (space separated)")

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

    return(record_dict)

# first, extract all info from the current lat long file and add to a dict
# then can check if new lat longs need to be added

current_lat_longs = {}

with open(args.current_lat_longs) as f:
    for line in f:
        line = line.strip()
        cols = line.split("\t")
        deme = cols[0]
        place = cols[1]
        lat_longs = (round(float(cols[2]), 1), round(float(cols[3]), 1))
        current_lat_longs[lat_longs] = {"deme": deme, "place": place}

# some bits of info in the genbank require some additional processing
# also will check if the lat longs already exist

def additional_processing(raw_record_dict):
    if raw_record_dict["country"] != "Unknown":
        country = raw_record_dict["country"].split(":")[0]
        state = raw_record_dict["country"].split(":")[1].split(",")[1].strip()
        town = raw_record_dict["country"].split(":")[1].split(",")[0].strip()
    else:
        state, town = "Unknown"

    raw_record_dict["country"] = country
    raw_record_dict["state"] = state
    raw_record_dict["town"] = town

    # also want to convert the date properly
    # useful info available here:
    # https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior
    date = raw_record_dict["collection_date"]
    python_date = datetime.strptime(date, "%b-%Y")
    formatted_date = str(python_date.year) + "-" + str(python_date.month) + "-XX"
    raw_record_dict["collection_date"] = formatted_date

    # check if the lat longs are in the current lat long file (rounded to 1 decimal place)
    lat = -round(float(raw_record_dict["lat_lon"].split(" ")[0]), 1)
    long = round(float(raw_record_dict["lat_lon"].split(" ")[2]), 1)

    if (lat, long) in current_lat_longs:
        print(lat, long)
    else:
        print("lat longs not found")
        print(lat, long)
        print(current_lat_longs)

    return(raw_record_dict)

# extract all required info from the genbank files

genbank_dict = {}

for genbank_file in args.genbank_files:
    original_dict = parse_genbank(genbank_file)
    genbank_dict[genbank_file] = additional_processing(original_dict)

print(genbank_dict)
