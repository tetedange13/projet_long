#/usr/bin/env python3

"""
Module gethering all functions linked both to the peeling process itself and
generation of the PU found
"""

import os
import subprocess as sub
import numpy as np
import multiprocessing as mp
import functools as ftls
import src.manage_io as mio
import src.external as ext


def peeling(peeled_pdb_path, peeled_pdb_id):
    """
    Run the process of peeling using the executable into bin/

    Args:
        peeled_pdb_path: Path (str) to pdb that has been peeled (or just the
        pdb to will be moved for alignment)
        peeled_pdb_id: Name (str) of the PDB that have been peeled

    Returns:
        A list of lines from the output of the peeling program, containing all
        the boundaries of the different PUs at each level of cutting
    """
    a_la_fac = False
    if a_la_fac:
        cmdLine_peel = ("bin/peel64 -pdb " + peeled_pdb_path +
                        " -dssp data/" + peeled_pdb_id + ".dss"
                        " -R2 95 -ss2 8 -lspu 20 -mspu 0 -d0 6.0 -delta 1.5"
                        " -oss 1 -p 0 -cp 0 -npu 16")

    else:
        cmdLine_peel = ("bin/peel32 " + peeled_pdb_path + " data/" +
                        peeled_pdb_id + ".dss 95 8 20 0 6.0 1.5 1 0 0")

    print("Peeling of " + peeled_pdb_id + " in progress...")
    outPeel_sub = sub.Popen(cmdLine_peel.split(),
                            stdout=sub.PIPE).communicate()[0]
    outPeel_w_dies = outPeel_sub.decode().split('\n')

    # Remove files generated and used by peeling execution:
    os.remove("data/" + peeled_pdb_id + '.dss')
    os.remove("file_ca_coo.pdb")
    os.remove("file_proba_contact.mat")
    if os.path.isfile("file_pu_delineation.mtx"):
        os.remove("file_pu_delineation.mtx")
        os.remove("file_matrix_pu_contact.mtx")

    list_out = [line for line in outPeel_w_dies if line and line[0] != '#']
    print("Peeling done!")
    
    return list_out


def peeled_to_dict(line):
    """
    Take a line from the peeling output and return a dictionary, where each key
    is the number of the PU (so starting at 1), with the value being the bounds
    of the considered PU

    Args:
        line: A line (str) from the output of the peeling program

    Returns:
        A dict with {PU_number:[inf_bound, sup_bound]}
    """
    all_bounds = list(map(int, line.split()[5:]))
    nb_tot_PU = len(all_bounds)/2
    i, nb_current_PU = 0, 1
    dict_PU = {}

    while nb_current_PU <= nb_tot_PU:
        dict_PU[nb_current_PU] = all_bounds[i:i+2]
        i += 2
        nb_current_PU += 1

    return dict_PU


def generate_PU_pdbs(dict_PU, level_cut, dict_coord_peeled, peeled_pdb_id):
    """
    Generate different pdb file, associated to each PU, based to the boundaries
    given as output of the peeling program

    Args:
        dict_PU: Dict containing the boundaries of each PU at a given level
        level_cut: The current level considered (int)
        dict_coord_peeled: Coordinates of the peeled pdb (dict)
        peeled_pdb_id: Name (str) of the PDB that have been peeled
    """
    nb_PU = len(dict_PU)

    for i in range(1, nb_PU+1):
        out_file = "results/" + peeled_pdb_id + "_PU_" + str(level_cut)
        with open(out_file + '_' + str(i) + '.pdb', 'w') as out_PU:
            inf_bound, sup_bound = dict_PU[i]

            for resID in range(inf_bound, sup_bound+1):
                out_PU.write(dict_coord_peeled[resID][1])


def get_bestAlgnd_PU(nb_PU, already_selcted, peeled_pdb_id, ref_pdb_id, level):
    """
    Align the different PU (that need to be aligned) against the reference pdb
    (using TMalign) and get the number of the PU that has the maximum TMscore

    Args:
        nb_PU: Total number of PU at this given level
        already_selcted: List of the PU that have already been aligned and
        chose as best aligned PU
        peeled_pdb_id: PDB ID (str) of the PDB that have been peeled
        ref_pdb_id: PDB ID (str) of the PDB to align against
        level: Current level (int) considered

    Returns:
        Index of the PU that is best aligned with the given PDB file
    """
    arr_scores = np.zeros(nb_PU, dtype=float) # To be able to use argmax()
    # Already aligned PU have their value set to -1 (will never be the max):
    for nb_PU_algnd in already_selcted:
        arr_scores[nb_PU_algnd-1] = -1

    for i in range(nb_PU):
        if (i+1) not in already_selcted:
            PU_name = peeled_pdb_id + "_PU_" + str(level) + '_' + str(i+1)
            # The "peel_longer" param is set to True, to avoid inversion
            # when the PUs are processed
            arr_scores[i] = ext.TM_align(PU_name, ref_pdb_id, True)

    return np.argmax(arr_scores) + 1


def write_algnd_PUs(PU_max_name, algnd_filename, idx):
    """
    Take the number of the best aligned PU, open the associated TM_file and
    append the file gethering all aligned PU with the coordinates of this
    new PU

    Args:
        PU_max_name: Filename (str) of the PU that had the best TMscore
        algnd_filename: Filename (str) of the file gethering all aligned PUs
        idx: The iteration (int) of PU alignment, to deal with opening mode of
        the "_PUs_algnd_" file
    """
    with open('results/' + PU_max_name + '.sup_all_atm', 'r') as sup_all_max:
        if idx == 0:
            aligned_PU = open('results/' + algnd_filename, 'w')
        else:
            aligned_PU = open('results/' + algnd_filename, 'a')

        for line in sup_all_max:
            if (line[0:4] == "ATOM") or ((line[0:6] == "HETATM") and
               ( (resName == "MET") or resName == "MSE") ):
               chain_ID = line[21:22].strip()

               if chain_ID == "A":
                   # Send chain A into the file getering all aligned PUs:
                   aligned_PU.write(line)

        aligned_PU.close()


def erase_algned(dict_coord_ref, ref_pdb_id, PU_max_name):
    """
    Write a new reference pdb file, by writing only the residues that have not
    been aligned yet

    Args:
        dict_coord_ref: Dict containing the coordinates of the reference pdb
        ref_pdb_id: Name (str) of the reference pdb to create
        PU_max_name: Filename (str) of the PU that had the best TMscore
    """
    # Get atoms of reference pdb that are aligned, they will be discarded
    # (corresponding to chain B of the .sup_atm file of the best aligned PU):
    with open('results/' + PU_max_name + '.sup_atm', 'r') as sup_max:
        set_to_discard = set()

        for line in sup_max:
            if (line[0:4] == "ATOM") or ((line[0:6] == "HETATM") and
               ( (resName == "MET") or resName == "MSE") ):
               chain_ID = line[21:22].strip()

               if chain_ID == "B":
                   set_to_discard.add(line[22:26].strip())

    # Now we write the new ref pdb, with aligned atoms erased:
    with open('results/' + ref_pdb_id + '.pdb', 'w') as pdb_file:
        resIDs = dict_coord_ref.keys()

        for resID in resIDs:
            resID_pdb, res_str = dict_coord_ref[resID]

            if resID_pdb not in set_to_discard:
                pdb_file.write(res_str)


def clean_sup_atm(peeled_pdb_id, level, nb_PU):
    """
    Clean .sup_atm and .sup_all_atm files after they have become useless

    Args:
        peeled_pdb_id: Name (str) of the PDB that have been peeled
        level: The current level considered (int)
        nb_PU: Total number of PU at this given level
    """
    for i in range(nb_PU):
            PU_name = peeled_pdb_id + "_PU_" + str(level) + '_' + str(i+1)

            if os.path.isfile("results/" + PU_name + ".sup_atm"):
                os.remove("results/" + PU_name + ".sup_atm")
            if os.path.isfile("results/" + PU_name + ".sup_all_atm"):
                os.remove("results/" + PU_name + ".sup_all_atm")


def get_best_level(res_peel, list_nb_PU):
    """
    Get the index of the level that is the best (trade-off between high TMscore
    and low number of PUs)

    Args:
        res_peel: Values (list) of the TMscores at each level
        list_nb_PU: List of the number of PUs associated with each level

    Returns:
        The index of the best level (i.e. with maximum ratio)
    """
    nb_levels = len(res_peel)

    if nb_levels == 1: # If only 1 level found for the protein peeled
        return 0

    else:
        ratio_score_nbPU = np.array(res_peel)/np.array(list_nb_PU)
        return np.argmax(ratio_score_nbPU)


def final_alignments(peeled_pdb_id, level, ref_pdb_id, peel_longer, out_peel):
    """
    Process the final step of alignment between all aligned PUs and the
    reference pdb (separated from the rest for checkpointing)

    Args:
        peeled_pdb_id: Name (str) of the PDB that have been peeled
        level: The current level considered (int)
        ref_pdb_id: Name (str) of the reference pdb to create
        peel_longer: Boolean telling if the peeled protein is longer (or not)
        than the reference protein
        out_peel: A list containing each line (str) from the peeling output

    Returns:
        The current level, the values of gdt-calculated TMscore, the (regular)
        TMscore and the total number of PU within the current level
    """
    PU_alignd_file = peeled_pdb_id + '-' + ref_pdb_id + '_PUs_' + str(level)
    nb_tot_PU = len(peeled_to_dict(out_peel[level-1]))

    TM_gdt = ext.gdt_pl(PU_alignd_file,
                        "results/" + ref_pdb_id + '_safe.pdb',
                        peel_longer)
    TMscore = ext.TM_score("results/" + ref_pdb_id + '_safe.pdb',
                           "results/" + PU_alignd_file + '.pdb',
                           peel_longer)

    return (level, TM_gdt, TMscore, nb_tot_PU)


def peel_calc(idx, out_peel, dictCoord_peeled, peeled_pdb_id,
            dictCoord_ref, ref_pdb_id, peel_longer):
    """
    Process to all the calculation and file manipulations for the
    "peeled-TMalignement"

    Args:
        idx: The number (int) of the iteration (used to deduce level)
        dictCoord_peeled: Dict containing the lines of the peeled pdb
        peeled_pdb_id: Name (str) of the PDB that have been peeled
        dictCoord_ref: Dict containing the lines of the reference pdb
        ref_pdb_id: Name (str) of the reference pdb to create
        peel_longer: Boolean telling if the peeled protein is longer (or not)
        than the reference protein

    Returns:
        The current level, the values of gdt-calculated TMscore, the (regular)
        TMscore and the total number of PU within the current level
    """
    level = idx + 1
    print("Proceeding peeling level", level)

    # Copy of the "safe" reference pdb, that will be proper to this level
    os.system("cp results/" + ref_pdb_id + "_safe.pdb results/" +
              ref_pdb_id + str(level) + '.pdb')

    dict_all_PU = peeled_to_dict(out_peel[idx])
    nb_tot_PU = len(dict_all_PU)

    generate_PU_pdbs(dict_all_PU, level, dictCoord_peeled, peeled_pdb_id)
    already_selcted = []

    # Then we loop on the number of PUs, to repeat the process
    for i in range(nb_tot_PU):
        nb_bestAlgnd_PU = get_bestAlgnd_PU(nb_tot_PU, already_selcted,
                                           peeled_pdb_id,
                                           ref_pdb_id + str(level),
                                           level)
        already_selcted.append(nb_bestAlgnd_PU)

        PUmax_name = (peeled_pdb_id + "_PU_" + str(level) + '_' +
                      str(nb_bestAlgnd_PU))
        algnd_filename = (peeled_pdb_id + '-' + ref_pdb_id + '_PUs_' +
                          str(level) + '.pdb')

        write_algnd_PUs(PUmax_name, algnd_filename, i)
        erase_algned(dictCoord_ref, ref_pdb_id + str(level), PUmax_name)
        clean_sup_atm(peeled_pdb_id, level, nb_tot_PU)

    # Remove pdb files of the PUs of the current level:
    for i in range(nb_tot_PU):
        os.remove("results/" + peeled_pdb_id + "_PU_" + str(level) + '_' +
                  str(i+1) + '.pdb')
    os.remove("results/" + ref_pdb_id + str(level) + '.pdb')

    return final_alignments(peeled_pdb_id, level, ref_pdb_id, peel_longer,
                            out_peel)


def peeled_TMalign(ref_pdb_path, ref_pdb_id, dictCoord_ref,
                   peeled_pdb_path, peeled_pdb_id, dictCoord_peeled,
                   peel_longer):
    """
    Function gethering the whole process of "peeled-TMalignment"

    Args:
        ref_pdb_path: Path (str) to the PDB to align against
        ref_pdb_id: Name (str) of the reference pdb to create
        dictCoord_ref: Dict containing the lines of the reference pdb
        peeled_pdb_path: Path (str) to pdb that has been peeled
        peeled_pdb_id: Name (str) of the PDB that have been peeled
        dictCoord_peeled: Dict containing the lines of the peeled pdb
        peel_longer: Boolean telling if the peeled protein is longer (or not)
        than the reference protein

    Returns:
        Three lists containing the different results
    """
    # We need a safe copy of the ref pdb, to reset it at each level:
    os.system("cp " + ref_pdb_path + " results/" + ref_pdb_id + "_safe.pdb")

    # Creation of dssp file (needed for peeling):
    if not os.path.isfile("data/" + peeled_pdb_id + ".dss"):
        os.system("bin/dssp64 -i " + peeled_pdb_path + " > data/" +
                  peeled_pdb_id + ".dss")

    # Peeling:
    out_peel = peeling(peeled_pdb_path, peeled_pdb_id)
    nb_tot_levels = len(out_peel)

    if True:
    # if not os.path.isfile("results/" + peeled_pdb_id + '_PUs_algnd_' +
    #                       str(nb_tot_levels) + '.pdb'):
        # Parallelized version (TO TIME):
        nb_cpu = mp.cpu_count() - 1
        # nb_cpu = 1
        my_pool = mp.Pool(nb_cpu)
        partial_func = ftls.partial(peel_calc,
                                      out_peel=out_peel,
                                      dictCoord_peeled=dictCoord_peeled,
                                      peeled_pdb_id=peeled_pdb_id,
                                      dictCoord_ref=dictCoord_ref,
                                      ref_pdb_id=ref_pdb_id,
                                      peel_longer=peel_longer)

        res_tot_peel = my_pool.map(partial_func, range(len(out_peel)))
        my_pool.close()

        # Serial version:
        # res_tot_peel = []
        # for idx in range(nb_tot_levels):
        #     res_tot_peel.append(peel_calc(idx, out_peel, dictCoord_peeled, peeled_pdb_id,
        #                          dictCoord_ref, ref_pdb_id, peel_longer))

    else:
        print("Found files of aligned PUs ! Skipping...")
        res_tot_peel = []
        for level in range(1, nb_tot_levels+1):
            res_tot_peel.append(final_alignments(peeled_pdb_id, level, ref_pdb_id,
                                          peel_longer, out_peel))

    # Putting outputs in proper lists:
    nb_levels = len(res_tot_peel)
    res_levels = [0.0] * nb_levels
    list_nb_PU = [0] * nb_levels
    res_gdt = [0.0] * nb_levels

    for res_tot in res_tot_peel:
        level, TM_gdt, TMscore, nb_tot_PU = res_tot
        res_levels[level-1] = TMscore
        list_nb_PU[level-1] = nb_tot_PU
        res_gdt[level-1] = TM_gdt

    # Cleaning useless "safe" file:
    os.remove('results/' + ref_pdb_id + '_safe.pdb')

    print('\n')
    return (res_levels, list_nb_PU, res_gdt)
