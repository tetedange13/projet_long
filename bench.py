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
import re
import urllib.request as urlreq
import pandas as pd
import pyquery as pyq


class AppURLopener(urlreq.FancyURLopener):
    version = "Mozilla/5.0"


def get_url_dom(sid_dom):
    """
    Get the exact url of the SCOP (domain) PDB, by looking through the
    HTML page associated to the domain
    """
    opener = AppURLopener()
    url_to_open = "http://scop.berkeley.edu/sid=" + sid_dom

    with opener.open(url_to_open) as response:
        pq = pyq.PyQuery(response.read())

    list_hrefs = []
    for a in pq.find('a'):
        current_href = a.get('href')
        if 'pdbstyle' in current_href and sid_dom in current_href:
            list_hrefs.append(current_href)

    if len(list_hrefs) != 1:
        print("ERROR: Several urls possible !! ")
        sys.error(2)
    else:
        return list_hrefs[0]


def dl_pdb(url_dom, pdb_id):
    """
    Download the given SCOP (domain) PDB file from the webpage
    """
    good_url = re.sub(r'(output=html)', 'output=txt', url_dom)

    print("Dowloading the good domain of " + pdb_id + ".pdb from the SCOP " +
          "website...")
    urlreq.urlretrieve(good_url, "data/" + pdb_id + '.pdb')
    print("Download finished !\n")




# MAIN:
if __name__ == "__main__":
    # Get dataset table:
    RIPC_txt = pd.read_csv('data/RIPC_dataset.txt', sep='\t')

    # Delete old csv file:
    if os.path.isfile('toto.csv'):
        os.remove("toto.csv")

    for idx, row in RIPC_txt.iterrows():
        len_dom1, len_dom2 = row['Length1'],row['Length2']
        dom1_sid, dom2_sid = row['Domain1'], row['Domain2']
        pdb_id_dom1, pdb_id_dom2 = dom1_sid[1:5], dom2_sid[1:5]
        chainID_sid1, chainID_sid2 = dom1_sid[5], dom2_sid[5]


        if idx == 1:
            print("LEN_DOM1 (PEEL):", len_dom1)
            print("LEN_DOM2 (REF):", len_dom2)

            # # Deals with the case where chainID is not specified ('_'):
            # if chainID_sid1 == '_':
            #     chainID_dom1 = 'A'
            # else:
            #     chainID_dom1 = chainID_sid1.upper()
            #
            # if chainID_sid2 == '_':
            #     chainID_dom2 = 'A'
            # else:
            #     chainID_dom2 = chainID_sid2.upper()

            # If the PDB is absent from the data/ folder, it is downloaded:
            if not os.path.isfile("data/" + pdb_id_dom1 + '.pdb'):
                url_dom1 = get_url_dom(dom1_sid)
                dl_pdb(url_dom1, pdb_id_dom1)
            if not os.path.isfile("data/" + pdb_id_dom2 + '.pdb'):
                url_dom2 = get_url_dom(dom2_sid)
                dl_pdb(url_dom2, pdb_id_dom2)

            print(pdb_id_dom1, pdb_id_dom2)
            os.system("./main.py -p data/" + pdb_id_dom1 + '.pdb -r data/' +
                      pdb_id_dom2 + ".pdb -b f")# -c " + chainID_dom1 + " -s " +
                      #chainID_dom2)

        # else:
        #     sys.exit()
