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


def TM_align(PU_name, pdbName_2):
    """
    Returns the TMscore?
    """
    # To view superimposed C-alpha traces of aligned regions: TM.sup
    #    To view superimposed C-alpha traces of all regions: TM.sup_all
    #    To view superimposed full-atom structures of aligned regions: TM.sup_atm
    #    To view superimposed full-atom structures of all regions: TM.sup_all_atm
    #    To view superimposed full-atom structures of all regions with ligands:
    #        TM.sup_all_atm_lig

    cmdLine_TM = ("bin/TMalign32 results/" + PU_name + '.pdb' + " results/" +
                   pdbName_2 + '.pdb' + " -o " + "results/" + PU_name + '.sup')
    out_TM = sub.Popen(cmdLine_TM.split(), stdout=sub.PIPE).communicate()[0]
    lines_TM = out_TM.decode()

    # Remove useless files:
    os.remove("results/" + PU_name + ".sup_all_atm_lig")
    os.remove("results/" + PU_name + ".sup_all")
    os.remove("results/" + PU_name + ".sup")
    #os.remove("results/" + PU_name + ".sup_all_atm")

    return float(lines_TM[1026:1033]) # The TMscore


def TM_to_pdb(pdb_name):
    """
    Open an output file from TMalign and extract the chain A, corresponding to the
    1st pdb now aligned on the reference pdb

    Args:
        pdb_name: Name (str) of the 1st pdb (the one to get the aligned coordinates)
    """
    with open("results/" + pdb_name + ".sup_all_atm") as all_atm_file, \
         open("results/" + pdb_name + '.pdb', 'w') as pdb_out:
        for line in all_atm_file:
            resName = line[17:20].strip()
            resID_pdb = line[22:26]

            if (line[0:4] == "ATOM") or ((line[0:6] == "HETATM") and
               ( (resName == "MET") or resName == "MSE") ):
               chain_ID = line[21:22].strip()
               if chain_ID == "A":
                   pdb_out.write(line)

    # Remove ".sup_all_atm" file (ow useless):
    os.remove("results/" + pdb_name + ".sup_all_atm")


def peeled_to_dict(line):
    """
    Take a line from the peeling output and return a dictionary, where each key
    is the number of the PU (so starting at 1), with the value being the bounds
    of the considered PU

    Args:
        line: A line (str) from thsplitted_line[5:]e output of peeling

    Returns:
        A dict with {PU_number:[inf_bound, sup_bound]}
    """
    splitted_line = [ int(elem) for elem in line.split()[4:] ]
    nb_tot_PU, all_bounds = splitted_line[0], splitted_line[1:]
    i, nb_current_PU = 0, 1
    dict_PU = {}

    while nb_current_PU <= nb_tot_PU:
        dict_PU[nb_current_PU] = all_bounds[i:i+2]
        i += 2
        nb_current_PU += 1

    return dict_PU


def generate_PU_pdbs(dict_PU, level_cut, dict_coord, name_pdb):
    """
    """
    nb_PU = len(dict_PU)

    for i in range(1, nb_PU+1):
        out_file = "results/" + name_pdb + "_PU_" + str(level_cut)
        with open(out_file + '_' + str(i) + '.pdb', 'w') as out_PU:
            inf_bound, sup_bound = dict_PU[i]

            for resID in range(inf_bound, sup_bound+1):
                out_PU.write(dict_coord[resID][1])


def generate_pdb_PU(bounds_PU, level_cut, dict_coord, name_pdb, idx):
    """
    """
    nb_PU = len(dict_PU)
    out_file =  name_pdb + "_PU_" + str(level_cut) + '_' + str(idx+1) + '.pdb'
    inf_bound, sup_bound = bounds_PU[(2*idx):(2*(idx+1))]

    with open("results/" + out_file, 'w') as out_PU:
        for resID in range(inf_bound, sup_bound+1):
            out_PU.write(dict_coord[resID][1])


def get_bestAlgnd_PU(dict_PU, already_selcted, pdbName_1, pdbName_2, level):
    """
    Aligned the different PU (that need to be aligned) against the

    Args:
        dict_PU: Dict containing the boundaries of each PU at a given level
        already_selcted: List of the PU that have already been aligned and
        chose as best aligned PU
        pdbName_1: Name (st                   r) of the PDB that have been peeled
        pdbName_2: Name (str) of the PDB to align against
        level: Current level (int) considered

    Returns:
        Index of the PU that is best aligned with the given PDB file
    """
    nb_PU = len(dict_PU)
    arr_scores = np.zeros(nb_PU, dtype=float) # To be able to use argmax()
    # Already aligned PU have their value set to -1 (will never be the max):
    for nb_PU_algnd in already_selcted:
        arr_scores[nb_PU_algnd-1] = -1

    for i in range(nb_PU):
        if (i+1) not in already_selcted:
            PU_name = pdbName_1 + "_PU_" + str(level) + '_' + str(i+1)
            arr_scores[i] = TM_align(PU_name, pdbName_2)

    print(arr_scores)
    nb_maxScore_PU = np.argmax(arr_scores) + 1

    # We can now delete all useless files:
    for i in range(nb_PU):
        if (i+1) != nb_maxScore_PU:
            PU_name = pdbName_1 + "_PU_" + str(level) + '_' + str(i+1)

            if os.path.isfile("results/" + PU_name + ".sup_atm"):
                os.remove("results/" + PU_name + ".sup_atm")
            if os.path.isfile("results/" + PU_name + ".sup_all_atm"):
                os.remove("results/" + PU_name + ".sup_all_atm")

    return nb_maxScore_PU


def test(bestAlgnd_PU, level):
    """
    Take the number of the best aligned PU, open the associated TM_file and with
    it:
        * Append the file gethering all aligned PU with the coordinats of this
        new PU
        * Get the resID of residues that need to be erased from the pdb (the one
        that PU are aligned against)

    Args:
        bestAlgnd_PU: The number (int) of the PU that has been determined as the
        one best aligned
        level: The current level considered (int)

    Returns:
        The ensemble (set) of residues that need to be erased (by their resID)

    """
    algnd_filename = 'PU_' + str(level) + '_algnd.pdb'
    with open('results/' + PU_name_max + '.sup_atm', 'r') as out_TM_max,\
         open('results/' + algnd_filename, 'a') as aligned_PU:
        set_to_discard = set()

        for line in out_TM_max:
            if (line[0:4] == "ATOM") or ((line[0:6] == "HETATM") and
               ( (resName == "MET") or resName == "MSE") ):
               chain_ID = line[21:22].strip()

               if chain_ID == "A":
                   # Send chain A into the file getering all aligned PUs:
                   aligned_PU.write(line)

               elif chain_ID == "B":
                   set_to_discard.add(line[22:26].strip())

        aligned_PU.write("TER\n")

    return set_to_discard


def erase_algned(dict_coord, set_to_discard, pdb_name):
    """
    Write a new reference pdb file, by writing only the residues that have not been
    aligned yet

    Args:
        dict_coord: Dict containing the coordinates of the reference pdb
        set_to_discard: The ensemble (set) of residues that have been aligned against
        the PU (so that must not to written)
        pdb_name: Name (str) of the reference pdb to create
    """
    with open('results/' + pdb_name + '.pdb', 'w') as pdb_file:
        resIDs = dict_coord.keys()

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
        # Peeling:idxMax+1
        cmdLine_peel = ("bin/peeling11_4.1 -pdb data/1aoh.pdb -dssp data/1aoh.dss"
                        " -R2 98 -ss2 8 -lspu 20 -mspu 0 -d0 6.0 -delta 1.5 "
                        "-oss 0 -p 0 -cp 0 -npu  16")
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

    # Creation of the results/ directory:
    if not os.path.isdir("results"):
        os.mkdir(results)

    # Copy pdb files towards results/ :
    os.system("cp data/" + pdbName_1 + ".pdb results/")
    os.system("cp data/" + pdbName_2 + ".pdb results/")

    # We start by making a simple TMalignment between both pdb:
    TM_align(pdbName_1, pdbName_2)
    TM_to_pdb(pdbName_1)


    level = 0
    for line in lines_peel:
        if line and line[0] != '#': # There is 1 element with empty str
            level += 1

            if level == 3:
                dict_all_PU = peeled_to_dict(line)
                generate_PU_pdbs(dict_all_PU, level, dictCoord_1, pdbName_1)

                already_selcted = []
                nb_tot_PU = len(dict_all_PU)

                # Then we loop on the number of PUs, to repeat the process
                for i in range(nb_tot_PU):
                #for i in range(2):
                    bestAlgnd_PU = get_bestAlgnd_PU(dict_all_PU, already_selcted,
                                                    pdbName_1, pdbName_2, level)
                    already_selcted.append(bestAlgnd_PU)
                    PU_name_max = pdbName_1 + "_PU_" + str(level) + '_' + str(bestAlgnd_PU)
                    print(PU_name_max)

                    set_to_discard = test(bestAlgnd_PU, level)

                    # Rewrite a pdb file (the reference one) with already aligned atoms
                    # deleted (corresponding to chain B):
                    erase_algned(dictCoord_2, set_to_discard, pdbName_2)

            #os.system()
