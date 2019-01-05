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
import pyquery as pyq


class AppURLopener(urlreq.FancyURLopener):
    version = "Mozilla/5.0"


def dl_pdb(pdb_id):
    """Download the given files from the webpage of a folder"""

    url_pdb = "https://files.rcsb.org/download/" + pdb_id + ".pdb"

    print("Dowloading " + pdb_id + ".pdb from the RCSB website...")
    urlreq.urlretrieve(url_pdb, "data/" + pdb_id + '.pdb')
    print("Download finished !\n")


def get_domain_range(sid_dom):
    opener = AppURLopener()

    my_url = "http://scop.berkeley.edu/sid=" + sid_dom

    with opener.open(my_url) as response:
        pq = pyq.PyQuery(response.read())

    titre = pq.find('title').text()
    print(titre)
    sys.exit()
    range_dom = titre.split(':')[-1]
    # return list(map(int, range_dom.spl))
    return range_dom




# MAIN:
if __name__ == "__main__":
    # Get dataset table:
    RIPC_txt = pd.read_csv('data/RIPC_dataset.txt', sep='\t')

    # Delete old csv file:
    if os.path.isfile('toto.csv'):
        os.remove("toto.csv")

    for idx, row in RIPC_txt.iterrows():
        dom1_sid, dom2_sid = row['Domain1'], row['Domain2']
        pdb_id_dom1, pdb_id_dom2 = dom1_sid[1:5], dom2_sid[1:5]
        chainID_sid1, chainID_sid2 = dom1_sid[5], dom2_sid[5]

        dom_range1 = get_domain_range(dom1_sid)
        dom_range2 = get_domain_range(dom2_sid)
        print(dom_range1, dom_range2)

        # if idx == 1:
        #     print("salut")
        #     # Deals with the case where chainID is not specified ('_'):
        #     if chainID_sid1 == '_':
        #         chainID_dom1 = 'A'
        #     else:
        #         chainID_dom1 = chainID_sid1.upper()
        #
        #     if chainID_sid2 == '_':
        #         chainID_dom2 = 'A'
        #     else:
        #         chainID_dom2 = chainID_sid2.upper()
        #
        #     # If the PDB is absent from the data/ folder, it is downloaded:
        #     if not os.path.isfile("data/" + pdb_id_dom1 + '.pdb'):
        #         dl_pdb(pdb_id_dom1)
        #     if not os.path.isfile("data/" + pdb_id_dom2 + '.pdb'):
        #         dl_pdb(pdb_id_dom2)
        #
        #     print(pdb_id_dom1+chainID_dom1, pdb_id_dom2+chainID_dom2)
        #     os.system("./main.py -p data/" + pdb_id_dom1 + '.pdb -r data/' +
        #               pdb_id_dom2 + ".pdb -c " + chainID_dom1 + " -s " +
        #               chainID_dom2 +" -b f")

        # else:
        #     sys.exit()
