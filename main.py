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


def generate_PU_pdbs(bounds_PU, level_cut, dict_coord, name_pdb):
    """
    """

    nb_PU = bounds_PU[0] # 1st element = total number of PU

    for i in range(nb_PU):
        out_file = "results/" + name_pdb + "_PU_" + str(level_cut) + '_' + str(i+1)
        with open(out_file + '.pdb', 'w') as out_PU:
            inf_bound, sup_bound = bounds_PU[(2*i+1):(2*(i+1)+1)]

            for resID in range(inf_bound, sup_bound+1):
                out_PU.write(dict_coord[resID])


def generate_pdb_PU(bounds_PU, level_cut, dict_coord, name_pdb, idx):
    """
    """

    out_file =  name_pdb + "_PU_" + str(level_cut) + '_' + str(idx+1) + '.pdb'
    inf_bound, sup_bound = bounds_PU[(2*idx):(2*(idx+1))]

    with open("results/" + out_file, 'w') as out_PU:
        for resID in range(inf_bound, sup_bound+1):
            out_PU.write(dict_coord[resID][1])


def TM_align(pdbName_1, pdbName_2, level, idx):
    """
    Returns the TMscore?
    """

    # To view superimposed C-alpha traces of aligned regions: TM.sup
    #    To view superimposed C-alpha traces of all regions: TM.sup_all
    #    To view superimposed full-atom structures of aligned regions: TM.sup_atm
    #    To view superimposed full-atom structures of all regions: TM.sup_all_atm
    #    To view superimposed full-atom structures of all regions with ligands:
    #        TM.sup_all_atm_lig


    PU_name = pdbName_1 + "_PU_" + str(level) + '_' + str(idx+1)
    cmdLine_TM = ("bin/32bits_TMalign results/" + PU_name + '.pdb' + " data/" +
                   pdbName_2 + '.pdb' + " -o " + "results/" + PU_name + '.sup')
    out_TM = sub.Popen(cmdLine_TM.split(), stdout=sub.PIPE).communicate()[0]
    lines_TM = out_TM.decode()
    #print(lines_TM)

    # Remove useless files:
    os.remove("results/" + PU_name + ".sup_all_atm_lig")
    os.remove("results/" + PU_name + ".sup_all")
    os.remove("results/" + PU_name + ".sup")
    os.remove("results/" + PU_name + ".sup_all_atm")

    return float(lines_TM[1026:1033]) # The TMscore


def erase_algned(dict_coord, set_to_discard, pdb_name):
    """
    """

    with open('results/' + pdb_name + '.pdb', 'w') as pdb_file:
        resIDs = dict_coord.keys()
        print(sorted([int(nb) for nb in set_to_discard]))
        for resID in resIDs:
            resID_pdb = dict_coord[resID][0]

            if resID_pdb not in set_to_discard:
                pdb_file.write(dict_coord[resID][1])



# MAIN:
if __name__ == "__main__":
    pdbName_1 = "1aoh"
    pdbName_2 = "1jlx"

    # Creation of dssp file:
    if not os.path.isfile("data/1aoh.dss"):
        os.system("bin/dssp data/1aoh.pdb > data/1aoh.dss")

    if not os.path.isfile("data/1jlx.dss"):
        os.system("bin/dssp data/1jlx.pdb > data/1jlx.dss")


    with open("data/" + pdbName_1 + ".pdb") as pdbFile_1, \
         open("data/" + pdbName_2 + ".pdb") as pdbFile_2:
        dictCoord_1 = mio.parse_pdb(pdbFile_1)
        dictCoord_2 = mio.parse_pdb(pdbFile_2)


    a_la_fac = False
    if a_la_fac:
        # Peeling:
        cmdLine_peel = ("bin/peeling11_4.1 -pdb data/1aoh.pdb -dssp data/1aoh.dss"
                        " -R2 98 -ss2 8 -lspu 20 -mspu 0 -d0 6.0 -delta 1.5 "
                        "-oss 0 -p 0 -cp 0 -npu 16")
        # The split function of the shlex module is used to generate a list of args:
        # os.system(cmd_line)
        #out, err = sub.Popen(shx.split(cmd_line), stdout=sub.PIPE).communicate()
        outPeel_1 = sub.Popen(cmdLine_peel.split(), stdout=sub.PIPE).communicate()[0]
        lines_peel = outPeel_1.decode().split('\n')

    else:
        with open("data/" + pdbName_1 + '_peeled.txt', 'r') as outPeel_1:
            content_peel = outPeel_1.read()
            print(content_peel)
            lines_peel = content_peel.split('\n')


    level = 0
    for line in lines_peel:
        #print(len(line))
        if line and line[0] != '#': # There is 1 element with empty str
            level += 1

            if level < 2:
                splitted_line = line.split()
                nb_PU = int(splitted_line[4])
                # bounds_PU = [int(limit) for limit in splitted_line[4:]]
                bounds_all_PU = [int(limit) for limit in splitted_line[5:]]
                arr_scores = np.zeros(nb_PU, dtype=float) # np array to use argmax()

                for i in range(nb_PU):
                    generate_pdb_PU(bounds_all_PU, level, dictCoord_1, pdbName_1, i)
                    arr_scores[i] = TM_align(pdbName_1, pdbName_2, level, i)

                idxMax = np.argmax(arr_scores)
                PU_name_max = pdbName_1 + "_PU_" + str(level) + '_' + str(idxMax+1)
                print(PU_name_max)
                algnd_filename = 'PU_' + str(level) + '_algnd.pdb'

                with open('results/' + PU_name_max + '.sup_atm', 'r') as out_TM_max,\
                     open('results/' + algnd_filename, 'a') as aligned_PU:
                    set_to_discard = set()

                    for line in out_TM_max:
                        if (line[0:4] == "ATOM") or ((line[0:6] == "HETATM") and
                           ( (resName == "MET") or resName == "MSE") ):
                           chain_ID = line[21:22].strip()

                           if chain_ID == "A":
                               aligned_PU.write(line)
                           elif chain_ID == "B":
                               set_to_discard.add(line[22:26].strip())

                    erase_algned(dictCoord_2, set_to_discard, pdbName_2)
                    # 1- Envoyer la chaine A (coord PU) dans le fichier rassemblant les coord
                    # des differentes PU alignees
                    # 2- Reecrire un fichier pdb (ref 1jlx) en suppr atoms associes a la B

    # From output of peeling:
    # bounds_all_PU = [1, 56, 57, 143, 144, 290]
    # nb_PU = 3
    # level = 2
