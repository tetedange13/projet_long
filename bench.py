#!/usr/bin/env python3

"""Benchmarking script
"""
# Usage:
#   bench.py -p <peelPdb> -r <refPdb> [-b <benchMode>]
#
# Options:
#   -h --help                  help
#   --version                  version of the script
#   -p --peelPdb = peeled_pdb  input pdb file, that will be peeled
#   -r --refPdb = ref_pdb      other input pdb file, that will be used as reference (not peeled)
#   -b --benchMode = benching  Mode benchmarking (bool) [default: False]
# """

import sys
import os
# import shutil
import urllib.request as urlreq
import pandas as pd


def dl_pdb(pdb_id):
    """Download the given files from the webpage of a folder"""

    url_pdb = "https://files.rcsb.org/download/" + pdb_id + ".pdb"

    print("Dowloading " + pdb_id + ".pdb from the RCSB website...")
    urlreq.urlretrieve(url_pdb, "data/" + pdb_id + '.pdb')
    print("Download finished !\n")


# MAIN:
if __name__ == "__main__":
    RIPC_txt = pd.read_csv('data/RIPC_dataset.txt', sep='\t')

    for idx, row in RIPC_txt.iterrows():
        dom1_sid, dom2_sid = row['Domain1'], row['Domain2']
        pdb_id_dom1, pdb_id_dom2 = dom1_sid[1:5], dom1_sid[1:5]
        chainID_sid1, chainID_sid2 = dom1_sid[5], dom2_sid[5]

        if idx < 5:
            # Deals with the case where chainID is not specified ('_'):
            if chainID_sid1 == '_':
                chainID_dom1 = 'A'
            else:
                chainID_dom1 = chainID_sid1.upper()

            if chainID_sid2 == '_':
                chainID_dom2 = 'A'
            else:
                chainID_dom2 = chainID_sid2.upper()

            # If the PDB is absent from the data/ folder, it is downloaded:
            if not os.path.isfile("data/" + pdb_id_dom1 + '.pdb'):
                dl_pdb(pdb_id_dom1)
            if not os.path.isfile("data/" + pdb_id_dom2 + '.pdb'):
                dl_pdb(pdb_id_dom2)

            # os.system("./main.py -p data/" + pdb_id_dom1 + '.pdb -r data/' +
            #           pdb_id_dom1 + " -c " + chainID_dom1 + " -s " +
            #           chainID_dom1 +" -b t")
            print("./main.py -p data/" + pdb_id_dom1 + '.pdb -r data/' +
                   pdb_id_dom1 + ".pdb -c " + chainID_dom1 + " -s " +
                 chainID_dom1 +" -b t")

        else:
            sys.exit()
