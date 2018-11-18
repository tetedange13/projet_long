#!/usr/bin/env python3


import numpy as np
import sys, os
import subprocess as sub
# import shlex as shx
import src.manage_io as mio


def calcul_distance(dico_pdb, dico_align):
    """
    Calcul distance between all residus pair by pair.
    Args: The dictionary containing the coordinate of the template's structure
    only CA and only the fragment which matched with the query
          The dictionary of one alignment
    Returns: A matrix of distance between all residus
    """
    size_query = len(dico_align['query'])
    # Create the empty matrix distance, which will contain distances
    # between all residus
    matrix_dist = np.empty( (size_query, size_query), dtype="float" )
    # Fill the matrix distance with NaN
    matrix_dist[:] = np.nan
    for i in range(0, size_query):
        # Only the half matrix is completed because it's symetric
        for j in range(i + 1, size_query):
            # A window of 5 residus is taken because it's not interesting to
            # calculate their distance, because they are too close.
            if j not in range(i, i + 5) :
                # If the two residus aren't a gap
                if (dico_align['query'][i] != "-" and
                    dico_align['query'][j] != "-" and
                    dico_align['template'][i] != "-"):
                    # calculate the distance between the two residus
                    dist = np.sqrt((dico_pdb[i]["x"] - dico_pdb[j]["x"])**2 +
                                   (dico_pdb[i]["y"] - dico_pdb[j]["y"])**2 +
                                   (dico_pdb[i]["z"] - dico_pdb[j]["z"])**2)
                    # fill the matrix with the calculated distance
                    matrix_dist[i, j] = dist
                else:
                    # -1 is added when one residu is a gap
                    matrix_dist[i, j] = -1
    return matrix_dist


def extract_subarrays(proba_arr, m): # PI(m)
    A = proba_arr[m:n, 0:m] # There is also np.ix_()
    B = proba_arr[0:m, m:n]
    C = proba_arr[0:m, 0:m]

    return (A, B, C) # Faut-il retourner plutot la somme des matrices ??
    # return (A*B - C**2)/((A+C)(B+C))


def calcul_energy(matrix_dist, dope_arr, index_aa, dico_align):
    """
    Calcul total energy of the query threaded on the template. A dope
    matrix is used to retrieve the energy by distance.
    Args: The matrix distance, which contains the distance between all residus
          The dope_arr is a 3D numpy array containing the associated energy by
    distance for two residus
          The index_aa is a dictionary containing the index of the dope_arr for
    each amino acid
          The dico_align is the dictionary of one alignment
    Returns: The total energy of the query threaded on the template
    """
    size_query = len(dico_align['query'])
    energy = 0
    for i in range(0, size_query):
        # Retrieve the first amino acid
        aa1 = dico_align['query'][i]
        # Only the half matrix is traveled
        for j in range(i + 1, size_query):
            # Retrieve the second amino acid
            aa2 = dico_align['query'][j]
            # Retrieve the distance calculated between the two residus
            dist = matrix_dist[i, j]
            # If the distance is -1, it's a gap so it's penalized
            if dist == -1:
                energy += cA.QUERY_GAP
            # No distance calculated because these two residus are too close
            # or they are too distant to have an impact
            elif not(np.isnan(dist) or dist > cA.CUTOFF):
                # Retrieve index of the amino acid in the dope matrix
                idx1 = index_aa[aa1]
                idx2 = index_aa[aa2]
                # Calculate the index in the list of distance in dope matrix
                # The distance begin to 0.25A and the step is 0.5
                idx_dist = int(round((dist - 0.25) / 0.5))
                # Retrieve the energy in the dope matrix between these two
                # residus with the distance calculated
                energy += dope_arr[idx_dist][idx1][idx2]
    return energy


def write_dict(dict_coord):
    pass


def generate_PU_pdbs(bounds_PU, level_cut, dict_coord, name_pdb):
    """
    """

    nb_PU = bounds_PU[0] # 1st element = total number of PU

    for i in range(nb_PU):
        out_file = "results/" + name_pdb + "_PU_" + str(level_cut) + '_' + str(i+1)
        with open(out_file + '.pdb', 'w') as out_PU:
            inf_bound, sup_bound = bounds_PU[(2*i+1):(2*(i+1)+1)]

            for resID in range(inf_bound, sup_bound+1):
                #for atoms in dict_coord[resID]:
                out_PU.write(dict_coord[resID])

def generate_pdb_PU(bounds_PU, level_cut, dict_coord, name_pdb, idx):
    """
    """

    out_file =  name_pdb + "_PU_" + str(level_cut) + '_' + str(idx+1) + '.pdb'
    inf_bound, sup_bound = bounds_PU[(2*idx):(2*(idx+1))]

    with open("results/" + out_file, 'w') as out_PU:
        for resID in range(inf_bound, sup_bound+1):
            out_PU.write(dict_coord[resID])


def TM_align(name_pdb, level, idx):
    """
    Returns the TMscore?
    """

    PU_filename = name_pdb + "_PU_" + str(level) + '_' + str(idx+1) + '.pdb'
    pdb_filename = name_pdb + '.pdb'
    cmdLine_TM = "bin/32bits_TMalign results/" + PU_filename + " data/" + pdb_filename
    out_TM = sub.Popen(cmdLine_TM.split(), stdout=sub.PIPE).communicate()[0]
    lines_TM = out_TM.decode()

    return float(lines_TM[1026:1033]) # The TMscore


# MAIN:
if __name__ == "__main__":
    # Creation of dssp file:
    if not os.path.isfile("data/1aoh.dss"):
        os.system("bin/dssp data/1aoh.pdb > data/1aoh.dss")

    a_la_fac = False
    if a_la_fac:
        # Peeling:
        cmdLine_peel = ("bin/peeling11_4.1 -pdb data/1aoh.pdb -dssp data/1aoh.dss"
                        " -R2 98 -ss2 8 -lspu 20 -mspu 0 -d0 6.0 -delta 1.5 "
                        "-oss 0 -p 0 -cp 0 -npu 16")
        # The split function of the shlex module is used to generate a list of args:
        # os.system(cmd_line)
        #out, err = sub.Popen(shx.split(cmd_line), stdout=sub.PIPE).communicate()
        out_peel = sub.Popen(cmdLine_peel.split(), stdout=sub.PIPE).communicate()[0]
        lines_peel = out_peel.decode().split('\n')

        pdb_name = "1aoh"
        with open("data/" + pdb_name + ".pdb") as pdb_file:
            dict_coord = mio.parse_pdb(pdb_file)

        level = 0
        for line in lines:
            #print(len(line))
            if line and line[0] != '#': # There is 1 element with empty str
                level += 1

                splitted_line = line.split()
                nb_PU = int(splitted_line[4])
                # bounds_PU = [int(limit) for limit in splitted_line[4:]]
                bounds_all_PU = [int(limit) for limit in splitted_line[5:]]

                for i in range(nb_PU):
                    bounds_all_PU[]
                    generate_pdb_PU(bounds_all_PU, level, dict_coord, name_pdb, i)
                    TM_align(pdb_name, level, i)

                #generate_PU_pdbs()

    with open("data/1aoh.pdb") as pdb_file:
        dict_coord = mio.parse_pdb(pdb_file)

    #bounds_PU = [2, 1, 143, 144, 290]
    bounds_all_PU = [3, 1, 56, 57, 143, 144, 290]
    generate_PU_pdbs(bounds_all_PU, 2, dict_coord, "toto")
    cmdLine_TM = "bin/32bits_TMalign results/toto_PU_2_1.pdb data/1aoh.pdb"
    out_TM = sub.Popen(cmdLine_TM.split(), stdout=sub.PIPE).communicate()[0]
    lines_TM = out_TM.decode()
    splitted_linesTM = lines_TM.split('\n')
    # salut = lines_TM.find("0.39161 (if normalized by length of Chain_2)")
    #lines_TM[1016:1033]
    TM_score = float(lines_TM[1026:1033])
    print(TM_score)
    # print(salut)
