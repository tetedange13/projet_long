#!/usr/bin/env python3

"""
This module deals with input and output files to the program.
It contains parsing and writing functions.
"""


def dope_to_arr(input_dope_file):
    """
    Args: The parsed dope_CAonly.par file, with only the potentials beetween
    the CA
    Returns: A 3D numpy array
    """
    nb_nrgies = 30
    nb_aa = 20
    dope_arr = np.zeros( (nb_nrgies, nb_aa, nb_aa), dtype="float" ) # 3D array
    i, j, k = 0, 0, 0
    index_aa = {} # To make correspond amino acids with index inside the array
    for one_line in input_dope_file:
        splitted_line = one_line.split()
        if k < 20: index_aa[three_to_one(splitted_line[2])] = k ; k +=1
        nrgies = [nrgy for nrgy in map(float, splitted_line[4:])]
        dope_arr[:, i, j] = nrgies
        if j != (nb_aa - 1): j += 1 # Next column
        else: i += 1 ; j = 0 # Next line and go back at the beginning of line
    return (dope_arr, index_aa)
    

def parse_atm_file(name_pdb, start):
    """
    Read a pdb or atm file
    Args: The name of the pdb or atm file
    Returns: A list of dictionaries containig coordinate of residus
    """
    list_dict = []
    res = []
    # Id for residu because in some pdb files the num residu
    # dosen't start to 1
    res_id = 1
    with open(name_pdb, "r") as filin:
        for line in filin:
            if line[0:6].strip() == "ATOM" and line[12:16].strip() == "CA":
                if res_id >= start:
                    list_dict.append({"residu":line[17:20].strip(),
                                 "numero":res_id,
                                 "x":float(line[30:38]),
                                 "y":float(line[38:46]),
                                 "z":float(line[46:54])})
                res_id += 1
    return list_dict

