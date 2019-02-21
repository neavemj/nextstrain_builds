#!/usr/bin/env python

# standardise headers for H9N2 sequences from NCBI
# these have already gone through Ivano's cleaning pipeline
# 21.02.2019

import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser("fix headers in H9N2 sequences downloaded from NCBI\n"
                                 "the sequences have already gone through Ivano's cleaning pipeline\n"
                                 "this script will make sure all metadata is present\n"
                                 "and also try and standardise the naming\n"
                                 "E.g., 'ck' => 'chicken'\n"
                                 "E.g., 'gadwall duck' => 'duck'\n")

parser.add_argument('-f', '--fasta_file', type = str,
        help = "fasta files to convert to nextstrain headers")
parser.add_argument('-o', '--output', type = str,
        help = "a name for the converted fasta file")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check required arguments are provided

if args.fasta_file is None:
    print("\n** required input missing\n"
          "** a fasta file is required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# the headers look like this
# >KJ000702 A/chicken/Shandong/07/2009 2009/07/05 4 (HA)
# However, if it is from human no host is listed (thus we have fewer columns)
# they all contain an NCBI accession and the initial A/
# I'll first work on the host bit
# grep ">" FilteredH9HAaligned_FWv1.fasta | cut -f 2- -d " " | cut -f 2 -d "/" | sort | uniq | wc
# 105
# Hmm, this is a lot - even I fix a few, Nextstrain could never dispaly so many colours
# might have have to go with just human or avian, or something
# some of the places seem like obvious spelling errors which I'll fix here

place_convert = {"FuJian": "Fujian", "GuangXi": "Guangxi", "GuiZhou": "Guizhou", "Heiilongjian": "Heilongjiang",
                 "HongKong": "Hong Kong", "Shaanxi": "Shanxi", "ShanXi": "Shanxi", "Shangdong": "Shandong",
                 "Yunnanbn": "Yunnan", "Yunnandh": "Yunnan", "Yunnanws": "Yunnan"}

with open(args.fasta_file) as fl:
    for record in SeqIO.parse(fl, "fasta"):
        print(record.id)
