#!/usr/bin/env python3

"""
This module deals with input and output files to the program.
It contains parsing and writing functions.
"""

# from Bio.PDB.Polypeptide import three_to_one, one_to_three


def parse_pdb_2(pdb_file):
    """
    Read a pdb or atm file
    Args: The name of the pdb or atm file
    Returns: A list of dictionaries containig coordinate of residus
    """

    dict_coord = {}
    # Id for residu because in some pdb files the num residu
    # dosen't start to 1
    resID = 1

    for line in pdb_file:

        if (line[0:4] == "ATOM") or ((line[0:6] == "HETATM") and
        ( (resName == "MET") or resName == "MSE") ):
        #if line[0:6].strip() == "ATOM" or "HETATM" and not "":
            resName = line[17:20].strip()
            atom_name = line[12:16].strip()
            dict_coord[resID] = {"resName":resName,
                                 "atom_name":atom_name,
                                 "num":line[22:26].strip(),
                                 "x":float(line[30:38]),
                                 "y":float(line[38:46]),
                                 "z":float(line[46:54])}

            if atom_name == "CA":
                resID += 1

    return dict_coord


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
