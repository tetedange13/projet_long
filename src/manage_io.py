#!/usr/bin/env python3

"""
This module deals with input and output files to the program.
It contains parsing and writing functions.
"""

import os


def extract_chain(pdb_path, chain_id='first'):
    """
    If the parameter "chain_id" is not specified, the default behaviour is todo
    extract the 1st chain from the pdb

    Args:
        pdb_path: Path (str) to the pdb file, from which the chain will be
        extracted
        chain_id: The ID (str, 1 or 2 letters) of the chain to extract

    Returns:
        The name of the file (with .pdb extension)
        The pdb id (without .pdb extension)
    """
    pdb_name = os.path.basename(pdb_path) # With .pdb extension
    pdb_id = os.path.splitext(pdb_name)[0] # Without .pdb extension
    out_filename = "results/" + pdb_id + chain_id + '.pdb'

    with open(pdb_path, 'r') as pdb_in, \
         open(out_filename, 'w') as pdb_out:
        i = 0

        for line in pdb_in:
            resName = line[17:20].strip()
            resID_pdb = line[22:26]

            if (line[0:4] == "ATOM") or ((line[0:6] == "HETATM") and
                        ( (resName == "MET") or resName == "MSE") ):
               chain_ID = line[21:22].strip()
               i += 1

               if chain_id == 'first':
                   if i == 1: # If it is the 1st ATOM line read
                       first_chain_id = chain_ID

                   if chain_ID == first_chain_id:
                       pdb_out.write(line)
                   else:
                       break

               else:
                   if chain_ID == chain_id:
                       pdb_out.write(line)

    if not os.stat(out_filename).st_size: # If the file is empty
        print("ERROR! The chain ID you specified does not belong to " +
              pdb_path + ' !\n')
        os.remove(out_filename) # Clean the empty file
        sys.exit(1)

    if chain_id == 'first':
        os.system("mv " + out_filename + " results/" + pdb_id +
                  first_chain_id + '.pdb')
        return (pdb_id + first_chain_id + '.pdb', pdb_id + first_chain_id)

    return (pdb_id + chain_id + '.pdb', pdb_id + chain_id)


def parse_pdb(pdb_file):
    """
    Read a pdb or atm file
    Args: The name of the pdb or atm file
    Returns: A list of dictionaries containig coordinate of residus
    """

    dict_coord = {}
    # Id for residu because in some pdb files the num residu
    # dosen't start to 1
    resID = 0

    for line in pdb_file:
        resName = line[17:20].strip()
        resID_pdb = line[22:26]

        if (line[0:4] == "ATOM") or ((line[0:6] == "HETATM") and
           ( (resName == "MET") or resName == "MSE") ):
            if line[12:16].strip() == "N": # Suppose that 1st = "N"
                resID += 1
                resID_pdb = line[22:26].strip() # Needed for erasing

            if resID not in dict_coord.keys():
                dict_coord[resID] = [resID_pdb, line]
            else:
                dict_coord[resID][1] += line

    return dict_coord
