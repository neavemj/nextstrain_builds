#!/usr/bin/env python

# copied this function from /augur/utils.py
# for some reason it doesn't split my meta fields correctly - want to find out why

import os, sys
import pandas as pd


def read_metadata(fname):
    if os.path.isfile(fname):
        metadata = pd.read_csv(fname, sep='\t' if fname[-3:]=='tsv' else ',',
                             skipinitialspace=True).fillna('')
        meta_dict = {}
        for ii, val in metadata.iterrows():
            if hasattr(val, "strain"):
                meta_dict[val.strain] = val.to_dict()
            elif hasattr(val, "name"):
                meta_dict[val.name] = val.to_dict()
            else:
                print("ERROR: meta data file needs 'name' or 'strain' column")

        return meta_dict, list(metadata.columns)
    else:
        print("ERROR: file with states does not exist")
        return {}, []


meta_return = read_metadata(sys.argv[1])

for record in meta_return[0]:
    print(meta_return[0][record]['platform'])
