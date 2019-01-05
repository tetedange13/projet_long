#!/usr/bin/env python3

"""
This module deals with input and output files to the program.
It contains parsing and writting functions.
"""

import os


def extract_chain(pdb_path, chain_id_arg='first'):
    """
    If the parameter "chain_id_arg" is not specified, the default behaviour is
    to extract the 1st chain from the pdb

    Args:
        pdb_path: Path (str) to the pdb file, from which the chain will be
        extracted
        chain_id_arg: The ID (str, 1 or 2 letters) of the chain to extract

    Returns:
        The name of the file (with .pdb extension)
        The pdb id (without .pdb extension)
    """
    pdb_name = os.path.basename(pdb_path) # With .pdb extension
    pdb_id = os.path.splitext(pdb_name)[0] # Without .pdb extension
    out_filename = "results/" + pdb_id + chain_id_arg + '.pdb'

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

               if chain_id_arg == 'first':
                   if i == 1: # If it is the 1st ATOM line read
                       first_chain_id = chain_ID

                   if chain_ID == first_chain_id:
                       pdb_out.write(line)
                   else:
                       break

               else:
                   if chain_ID == chain_id_arg:
                       pdb_out.write(line)

    if not os.stat(out_filename).st_size: # If the file is empty
        print("ERROR! The chain ID you specified does not belong to " +
              pdb_path + ' !\n')
        os.remove(out_filename) # Clean the empty file
        sys.exit(1)

    # Rename file, in the case where no chain had been specified:
    if chain_id_arg == 'first':
        os.system("mv " + out_filename + " results/" + pdb_id +
                  first_chain_id + '.pdb')
        return (pdb_id + first_chain_id + '.pdb', pdb_id + first_chain_id)

    return (pdb_id + chain_id_arg + '.pdb', pdb_id + chain_id_arg)


def parse_pdb(pdb_file):
    """
    Read a pdb and get all the lines associated with each residue
    The keys are the residue ID, but reindexed and the value contains the old
    residue ID (the original one, from the PDB)

    Args:
        pdb_file: The name (str) of the pdb or atm file

    Returns:
        A dict where values are all the lines associates to a given residue and
        the size of the read PDB
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

    # The resID is now equivalent to the size of the protein given as argument
    return (resID, dict_coord)


def write_to():
    pass
