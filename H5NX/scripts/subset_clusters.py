#!/usr/bin/env python

# take clstr file from CD-HIT, fasta file and meta file
# then subset according to wanted criteria
# e.g. discard seqs > 99%, same town, same host, etc.
# 29.04.2019

import sys
import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser("subset sequences based on clustering and meta data\n")

parser.add_argument('-f', '--fasta', type = str,
        help = "fasta file with sequences to subset")
parser.add_argument('-m', '--meta', type = str,
        help = "meta file for sequences to subset")
parser.add_argument('-c', '--clusters', type = str,
        help = ".clstr file from a CD-HIT run")
parser.add_argument('-o', '--fasta_output', type = str,
        help = "a name for the subsetted fasta file")
parser.add_argument('-p', '--meta_output', type = str,
        help = "a name for the subsetted meta file")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check required arguments are provided

if args.fasta is None or \
   args.meta is None or \
   args.clusters is None or \
   args.fasta_output is None or \
   args.meta_output is None:
    print("\n** required input file is missing\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# create a dictionary containing each strain and its meta data
# I'll keep the whole line as an entry for easy writing out later

meta_dict = {}

with open(args.meta) as fl:
    meta_header = next(fl)
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        strain = cols[0]
        meta_dict[strain] = {
        "line": line,
        "location": cols[5],
        "country": cols[6],
        "date": cols[7],
        "host": cols[8]
        }

# use the clustering results to subset similar sequences
# will create two dictionaries for easier searching

cluster_to_strain = {}
strain_to_cluster = {}

with open(args.clusters) as fl:
    for line in fl:
        line = line.strip()
        if line.startswith(">"):
            cluster = line.lstrip(">").replace(" ", "_")
        else:
            cols = line.split(",")
            strain = cols[1].strip().lstrip(">").split("...")[0].strip()

            if cluster in cluster_to_strain:
                cluster_to_strain[cluster].add(strain)
            else:
                cluster_to_strain[cluster] = {strain}

            strain_to_cluster[strain] = cluster


# now go through fasta and check for duplicates

total_strains = 0
records_to_write = []
record_ids_to_write = []
strains_to_write = set()
meta_to_write = {}
duplicate_ids = 0
duplicate_id_list = []
singleton_clusters = 0
multiple_clusters = 0
strains_in_multiple_clusters_written = 0
dropped_strains = 0
processed_clusters = set()

with open(args.fasta) as fl:
    for record in SeqIO.parse(fl, "fasta"):
        total_strains += 1
        record_id = record.id
        record_cluster = strain_to_cluster[record_id]
        # records that are the single member of a cluster should be written
        # these are not similar to any other strain
        if len(cluster_to_strain[record_cluster]) == 1:
            singleton_clusters += 1
            if record_id not in record_ids_to_write:
                records_to_write.append(record)
                record_ids_to_write.append(record_id)
                meta_to_write[record_id] = meta_dict[record_id]["line"]
            else:
                duplicate_ids += 1
                duplicate_id_list.append(record_id)
        else:
            # if a cluster has already been processed, several strains might need to be
            # written independantly of the below code
            if record_id in strains_to_write:
                records_to_write.append(record)
                record_ids_to_write.append(record_id)
                meta_to_write[record_id] = meta_dict[record_id]["line"]
            # if the strains in a cluster have alread been checked, skip remaining code
            if record_cluster in processed_clusters:
                continue
            else:
                # this first strain in a cluster will become the representative
                # other stains in the cluster will be compared only to this
                # could change this to be a bit smarter
                processed_clusters.add(record_cluster)
                records_to_write.append(record)
                record_ids_to_write.append(record_id)
                meta_to_write[record_id] = meta_dict[record_id]["line"]
                strains_in_multiple_clusters_written += 1
                # i'll make sets of the meta data fields
                # then I can add to these as each isolate goes through the loop
                # this ensures that isolates are checked against all e.g. countries present
                location_set = {meta_dict[record_id]["location"]}
                country_set = {meta_dict[record_id]["country"]}
                date_set = {meta_dict[record_id]["date"]}
                host_set = {meta_dict[record_id]["host"]}
                # now check if we want to write any other records in this cluster
                for strain in cluster_to_strain[record_cluster]:
                    multiple_clusters += 1
                    # dont loop through the representaive strain above
                    if strain == record_id: continue
                    # use representative strain info to check against other strains
                    # if any of these fields are different, the stain will be kept
                    if any([ \
                       meta_dict[strain]["location"] not in location_set, \
                       meta_dict[strain]["country"] not in country_set, \
                       meta_dict[strain]["host"] not in host_set, \
                       meta_dict[strain]["date"] not in date_set \
                       ]):

                       # one of these must be different to the current cluster metadata
                       # add the new entries
                       location_set.add(meta_dict[strain]["location"])
                       country_set.add(meta_dict[strain]["country"])
                       host_set.add(meta_dict[strain]["host"])
                       date_set.add(meta_dict[strain]["date"])

                       # can't add strain directly to my record list
                       # because not looping through that yet
                       # will add the id to a list instead, then check against this as records come through
                       strains_in_multiple_clusters_written += 1
                       strains_to_write.add(strain)
                    else:
                        dropped_strains += 1


# write list of duplicate ids, and fasta and meta files

with open("duplicate_ids.txt", "w") as fl:
    for dup in duplicate_id_list: fl.write(dup + "\n")

SeqIO.write(records_to_write, args.fasta_output, "fasta")

with open(args.meta_output, "w") as fl:
    fl.write(meta_header)
    for strain_id in meta_to_write:
        fl.write(meta_to_write[strain_id] + "\n")


print("There were {} strains in the original fasta file".format(total_strains))
print("There were {} clusters formed".format(len(cluster_to_strain)) + "\n")

print("There were {} isolates that were 'singletons' and did not cluster with any other sequence".format(singleton_clusters))

print("There were {} isolates that clustered with others".format(multiple_clusters))
print("There were {} isolates that clustered with others and were written".format(strains_in_multiple_clusters_written))

print("There were {} isolates that were dropped due to the similarity filtering".format(dropped_strains))
print("There were {} isolates that were dropped due to duplicate strain ids".format(duplicate_ids))
print("\t* the ids for these were written to the file 'duplicate_ids.txt'" + "\n")

print("\nThere were {} total records written".format(len(records_to_write)))
print("There were {} total meta data records written".format(len(meta_to_write)))








