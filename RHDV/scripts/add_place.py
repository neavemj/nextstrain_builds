#!/usr/bin/env python

# need to summarise / collapse some of the gps co-ordinates
# Nextstrain won't accept more than 180 places for the geo bit
# Matthew J. Neave 12.12.2018

# will round co-ordinates to 1 decimal place

import sys

coords_file = sys.argv[1]
new_coords_file = open(sys.argv[2], "w")
meta_files = sys.argv[3:]

# first make dictionary of the exact GPS co-ordinates

exact_coords = {}

with open(coords_file) as f:
    next(f)
    for line in f:
        line = line.strip()
        cols = line.split("\t")
        strain = cols[0].strip()
        lat = float(cols[1].strip())
        lat_rounded = round(lat, 1)
        long = float(cols[2].strip())
        long_rounded = round(long, 1)
        exact_coords[strain] = {"lat": lat, "lat_rounded": lat_rounded, "long": long, "long_rounded": long_rounded}

# now create places from new 'rounded' coords_file

places = 0
lat_long_place = {}
new_meta = {}

for strain in exact_coords:
    new_coords = str(exact_coords[strain]["lat_rounded"]) + " " + str(exact_coords[strain]["long_rounded"])
    if new_coords in lat_long_place:
        new_meta[strain] = lat_long_place[new_coords]
    else:
        places += 1
        place = "place_" + str(places)
        lat_long_place[new_coords] = place
        new_meta[strain] = lat_long_place[new_coords]

# write new coords file

for coord in lat_long_place:
    new_coords_file.write("\t".join(["place", lat_long_place[coord], coord.split(" ")[0], coord.split(" ")[1]]) + "\n")

# now go through meta data files and add 'place'

for meta_fl in meta_files:
    with open(meta_fl) as fl:
        raw_header = next(fl)
        new_header = raw_header.strip() + "\tplace\n"
        new_meta_fl = open(meta_fl.rstrip(".meta.tsv") + ".meta.place.tsv", "w")
        new_meta_fl.write(new_header)
        for line in fl:
            line = line.strip()
            cols = line.split()
            strain = cols[0].strip()
            try:
                place = new_meta[strain]
            except:
                print("\nWarning: having trouble processing the following meta-data line:")
                print(line + "\n")
            new_meta_fl.write(line + "\t" + place + "\n")
