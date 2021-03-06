#!/usr/bin/env python

# standardise headers for H9N2 sequences from NCBI
# these have already gone through Ivano's cleaning pipeline
# 21.02.2019

import sys
import argparse
import re
from Bio import SeqIO


parser = argparse.ArgumentParser("fix headers in H9N2 sequences downloaded from NCBI\n"
                                 "the sequences have already gone through Ivano's cleaning pipeline\n"
                                 "this script will make sure all metadata is present\n"
                                 "and also try and standardise the naming\n"
                                 "'gadwall duck' => 'avian'\n"
                                 "'water' => 'environmental'\n"
                                 "\n*strains will be dropped if they do not meet nextstrain requirements\n")

parser.add_argument('-f', '--fasta_file', type = str,
        help = "fasta files to convert to nextstrain headers")
parser.add_argument('-g', '--geo', type = str,
        help = "geo_synonyms.tsv file to convert label to country, division, location")
parser.add_argument('-o', '--fasta_output', type = str,
        help = "a name for the converted fasta file")
parser.add_argument('-m', '--meta_output', type = str,
        help = "a name for the meta file")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check required arguments are provided

if args.fasta_file is None or \
   args.geo is None or \
   args.fasta_output is None or \
   args.meta_output is None:
    print("\n** required input file is missing\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# the headers look like this
# >KJ000702 A/chicken/Shandong/07/2009 2009/07/05 4 (HA)
# However, if it is from human no host is listed (thus we have fewer columns)
# they all contain an NCBI accession and the initial A/

def correct_strain_format(strain):
    # regular expression to get the strain name bit
    # must start with A, then any letters, then letters or numbers repeating
    # can contain dashes, underscores, spaces, dots, apostrophies, or brackets
    # then end with a series of digits (year). Can be 2009 or 09
    # note the space at the end of the match must be there
    #match = re.match(r'A/[A-Za-z-_ \.\']+/([A-Za-z0-9-_ \.\'\(\)]+/){0,3}[0-9]+ ', strain)
    match = re.match(r'A/[A-Za-z-_ \.\']+/([A-Za-z0-9-_ \.\'\(\)]+/){0,3}[0-9]+ ', strain)
    if match:
        return match.group(0)
    else:
        print("Filtered due to incorrect strain format:", strain)
        return False

# lets check the host bit
# grep ">" FilteredH9HAaligned_FWv1.fasta | cut -f 2- -d " " | cut -f 2 -d "/" | sort | uniq | wc
# 105
# Hmm, this is a lot - even I fix a few, Nextstrain could never display so many colours
# might have have to go with just human, avian, or environmental - this is what the nextstrain guys do

avian_host = ['ruddy turnstone', 'rosy-billed pochard', 'common teal', 'knot', 'turtledove', 'shorebird', 'muscovy duck', 'silkie chicken', 'green winged teal', 'wild waterfowl', 'pigeon', 'greater white-fronted goose', 'poultry', 'crow', 'coot', 'falcon', 'sparrow', 'gull', 'cackling goose', 'pelican', 'Chicken', 'goose', 'sanderling', 'Muscovy duck', 'bewicks swan', 'chiken', 'partridge', 'teal', 'Chinese francolin', 'avian', 'brambling', 'blue-winged teal', 'northern pintail', 'Chinese Hwamei', 'baikal teal', 'red knot', 'quail', 'parakeet', 'wild chicken', 'ck', 'Anas platyrhynchos', 'duck', 'Korean native chicken', 'black chicken', 'emperor goose', 'Anser fabalis', 'stone curlew', 'common murre', 'Himantopus himantopus', 'wild bird', 'turkey', 'mallard duck', 'chicken', 'green-winged teal', 'gadwall', 'northern shoveler', 'bird', 'laughing gull', 'Duck', 'mallard', 'common redshank', 'black-headed gull', 'egret', 'Mediterranean gull', 'bean goose', 'guineafowl', 'American oystercatcher', 'black-billed magpie', 'Avian', 'gadwall duck', "Bewick's swan", 'African Stonechat', 'ostrich', 'Japanese Quail', 'pink-footed goose', 'pheasant', 'American wigeon', 'Eurasian wigeon', 'broiler chicken', 'dove', 'Anthropoides virgo', 'American green-winged teal', 'white-fronted goose']

environment_host = ['feces', 'environment', 'wild bird feces', 'pigeon feces']

nonhuman_mammal_host = ['canine', 'swine', 'pika', 'equine', 'mink', 'Sw', 'swine']

other_host = ['ferret', 'insect']

def format_host(host):
    '''
    Fix host formatting
    Change any bird species to avian
    '''
    if host is not None:

        if host in avian_host:
            host_type = "avian"
            return(host_type)

        if host in environment_host:
            host_type = "environment"
            return(host_type)

        if host in nonhuman_mammal_host:
            host_type = "nonhuman_mammal"
            return(host_type)

        if host in other_host:
            host_type = "other"
            return(host_type)

        else:
            print("Could not assign host {}".format(host))
            print("Leaving current host designation as is")
            return(host)


# some of the places seem like obvious spelling errors which I'll fix here
# YN might be a province in China - will leave it for the moment - it will get filtered in the next step

place_convert = {"FuJian": "Fujian", "GuangXi": "Guangxi", "GuiZhou": "Guizhou", "Heiilongjian": "Heilongjiang", "Heiilongjiang": "Heilongjiang",
                 "HongKong": "Hong Kong", "Shaanxi": "Shanxi", "ShanXi": "Shanxi", "Shangdong": "Shandong",
                 "Yunnanbn": "Yunnan", "Yunnandh": "Yunnan", "Yunnanws": "Yunnan", "DE": "Delaware", "Emirates": "United Arab Emirates",
                 "Heibei": "Hebei", "HeBei": "Hebei", "Interior Alaska": "Alaska", "MN": "Minnesota", "NZL": "New Zealand",
                 "PT": "Portugal", "ShanDong": "Shandong", "Southcentral Alaska": "Alaska", "TX": "Texas", "YN": "YN", "Korea": "South Korea",
                 "Delaware Bay": "Delaware"}

# in the nextstrain github, there is a big list of places and what country they are in
# they also have a decent list of lat longs
# i'll use these to check if I can find the places and country

country_set = set()
place_to_country = {}
place_to_country_camel = {}

with open(args.geo) as fl:
    for line in fl:
        line = line.strip()
        if line.startswith("#"):
            continue
        if line:
            cols = line.split("\t")
            # change these from camel case to underscore separated
            country = cols[1]
            country_unders = "_".join(re.findall('[A-Z][^A-Z]*', country))
            location = cols[3]
            location_unders = "_".join(re.findall('[A-Z][^A-Z]*', location))
            country_set.add(country_unders)
            place_to_country[location_unders] = country_unders
            place_to_country_camel[location] = country_unders

# now read through fasta file and check headers
# create metadata file for Nextstrain

filtered_incorrect_strain_format = 0
filtered_host_first_column = 0
filtered_no_country_latlongs = 0
records_to_write = []
meta_to_write = {}

with open(args.fasta_file) as fl:
    for record in SeqIO.parse(fl, "fasta"):
        desc = record.description
        ncbi_id = desc.split(" ")[0]
        strain_description = " ".join(desc.split(" ")[1:])
        strain_id = correct_strain_format(strain_description)
        if not strain_id:
            filtered_incorrect_strain_format += 1
            continue
        if strain_id:
            id_components = strain_id.split("/")

            # headers with 5 fields tend to more correct
            if len(id_components) == 5:
                host = id_components[1]
                host_type = format_host(host)
                place = id_components[2]
                year = id_components[4].strip()

            # had more problems with headers with only 4 fields
            elif len(id_components) == 4:
                # theoretically host should always be human in this case
                host = "human"
                host_type = "human"
                place = id_components[1]
                # this is unfortunately not always true
                # sometimes component[1] is the host
                if place in avian_host or \
                place in environment_host or \
                place in nonhuman_mammal_host or \
                place in other_host:
                    # this is only 4 sequences
                    # I'll filter these out for the moment
                    print("Filtered due to 4 fields, but host in first column:", strain_id)
                    filtered_host_first_column += 1
                    continue
                year = id_components[3].strip()

            # the new Indonesian isolates have a header with 6 fields
            elif len(id_components) == 6:
                host = id_components[1]
                host_type = format_host(host)
                place = id_components[2]
                year = id_components[5].strip()

            # ok, hopefully the host, place, and year fields have been retrieved
            # will do a bit of formatting on these fields

            # fix spelling mistakes or abbreviations in the place
            if place in place_convert:
                place = place_convert[place]

            # for some reason the places in the geo_synonyms have underscores instead of spaces
            place = place.replace(" ", "_")

            # try and assign a country to a place
            # in some cases, the place given might already be the country
            if place in country_set:
                country = place
            elif place in place_to_country:
                    country = place_to_country[place]
            elif place in place_to_country_camel:
                    country = place_to_country_camel[place]
            else:
                print("Filtered because could not assign to a country:", strain_id, place)
                filtered_no_country_latlongs += 1
                continue

            # now deal with the year of collection
            # sometimes this is 4 digits (1997), sometimes 2 digits (97)
            abbrev_year = re.match(r'[0-9][0-9]$', year)
            if abbrev_year:
                abbrev_year = abbrev_year.group(0)
                if int(year) < 40:
                    year = "20" + abbrev_year
                else:
                    year = "19" + abbrev_year
            date = year + "-XX-XX"

            # want to have the strain id as the fasta header (not NCBI)
            strain_id = strain_id.strip().replace(" ", "_")
            record.id = strain_id
            record.description = ""
            records_to_write.append(record)
            meta_to_write[strain_id] = [strain_id, "AIV", ncbi_id, "HA", "H9N2", place, country, date, host, host_type, "unpublished", "unpublished", "unpublished"]



# now write out fasta and meta files
# meta headers should be:
# strain virus accession segment type place country date host authors title journal
# A/chicken/Henan AIV M542245 HA H9N2 Henan China 2018 mallard avian unpublished unpublished unpublished
SeqIO.write(records_to_write, args.fasta_output, "fasta")

with open(args.meta_output, "w") as fl:
    fl.write("\t".join(["strain", "virus", "accession", "segment", "type", "location", "country", "date", "host", "host_type", "authors", "title", "journal"]) + "\n")
    for strain_id in meta_to_write:
        fl.write("\t".join(meta_to_write[strain_id]) + "\n")


print("\nTotal fasta records written {}\n"
      "Total meta data entries: {}\n"
      "The above numbers should theoretically match!!\n\n"
      "Filtered due to incorrect strain id: {}\n"
      "Filtered due to host in first column: {}\n"
      "Filtered due to missing country: {}\n".format(len(records_to_write), len(meta_to_write), filtered_incorrect_strain_format, filtered_host_first_column, filtered_no_country_latlongs))













