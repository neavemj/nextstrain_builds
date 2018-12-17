#!/usr/bin/env python

# Kim Blasdell sent me the coords for each town in their study
# want to add this info to the metadata files
# plus create a lat long file for Cooktown
# also need to create 'place' for existing lat long in other files

import sys

blas_file = sys.argv[1]
new_lat_longs = open(sys.argv[2], "w")
meta_files = sys.argv[3:]

blas_dict = {}

with open(blas_file) as f:
    next(f)
    for line in f:
        line = line.strip()
        cols = line.split("\t")
        town = cols[4].replace(" ", "_")
        isolate = cols[6]
        try:
            lat = cols[11]
            long = cols[12]
        except:
            print("lat long not found for: {}".format(isolate))
            continue
        blas_dict[isolate] = {"town": town, "lat": lat, "long": long}

# now add this new info to the meta files

place_lat_longs = {}
no_lat_longs = set()

for meta_file in meta_files:
    with open(meta_file) as f:
        new_name = meta_file.split(".")[0] + "_meta.new.tsv"
        new_meta_file = open(new_name, "w")
        header = next(f).strip()
        new_meta_file.write(header + "\tplace\n")
        for line in f:
            line = line.strip()
            cols = line.split("\t")
            strain = cols[0]
            if strain in blas_dict:
                new_meta_file.write(line + "\t" + blas_dict[strain]["town"] + "\n")
            else:
                no_lat_longs.add(strain)

# write out the lat long files

written_towns = []

for strain in blas_dict:
    if blas_dict[strain]["town"] in written_towns:
        continue
    else:
        written_towns.append(blas_dict[strain]["town"])
        new_lat_longs.write("place\t" + blas_dict[strain]["town"] + "\t" \
        + blas_dict[strain]["lat"]+ "\t" + blas_dict[strain]["long"] + "\n")
