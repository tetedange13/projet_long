#!/usr/bin/env python3

"""PeelAlign script
Usage:
  main.py -p <peelPdb> -r <refPdb>

Options:
  -h --help                  help
  --version                  version of the script
  -p --peelPdb = peeled_pdb  input pdb file, that will be peeled
  -r --refPdb = ref_pdb      other input pdb file, that will be used as reference (not peeled)
"""


from docopt import docopt
import numpy as np
import sys, os
import subprocess as sub
import src.manage_io as mio
import re
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import src.peeling as pl


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


def TM_score(peeled_pdb_path, ref_pdb_path):
    """
    Using the TMscore binary, calculate the TMscore of a given pdb, containing
    all aligned PUs and the reference pdb (so no alignment is processed)

    Args:
        peeled_pdb_path: Path (str) to pdb that has been peeled (or just the
        pdb to will be moved for alignment)
        ref_pdb_path: Path (str) to the PDB to align against

    Returns:
        The value of the associated TMscore
    """
    cmdLine_TM = ("bin/TMscore32 " + peeled_pdb_path + " " + ref_pdb_path)
    out_TM = sub.Popen(cmdLine_TM.split(), stdout=sub.PIPE).communicate()[0]
    lines_TM = out_TM.decode()
    # print(lines_TM)

    regex_TMscore = re.compile("(?:TM-score.+= )([0-9]\.[0-9]*)[ $]")
    searchObj = re.search(regex_TMscore, lines_TM)

    # It is possible to have a case where the TMscore does not find any
    # residues in common, so we return -1
    if searchObj:
        return float(searchObj.group(1))
    return -1


def parMATT(peeled_pdb_path, ref_pdb_path):
    cmdLine_parMatt = ("bin/parMATT/bin/parMatt32 " + peeled_pdb_path + " " +
                       ref_pdb_path + " -t 1 -o output")
    out_parMatt = sub.Popen(cmdLine_parMatt.split(),
                            stdout=sub.PIPE).communicate()[0]
    # print(out_parMatt.decode())

    PEELED_PDB, PEELED_PDB_ID = extract_chain("output.pdb", "A")
    REF_PDB, REF_PDB_ID = extract_chain("output.pdb", "B")

    # Remove useless files produced by parMATT:
    for ext in ('pdb', 'spt', 'fasta', 'txt'):
        os.remove("output." + ext)

    return TM_score("results/" + PEELED_PDB, "results/" + REF_PDB)


def TM_align(PU_name, ref_pdb_name):
    """
    Using the TMalign binary, proceed to a structural alignment between a given
    PU and the reference pdb (TMscore is returned as the stdout of TMalign
    program)

    Args:
        PU_name: Name (str) of the PU to align
        ref_pdb_name: Name (str) of the PDB to align against

    Returns:
        The value of the associated TMscore
    """
    cmdLine_TM = ("bin/TMalign32 results/" + PU_name + '.pdb' + " results/" +
                  ref_pdb_name + '.pdb' + " -o " + "results/" + PU_name + '.sup')
    out_TM = sub.Popen(cmdLine_TM.split(), stdout=sub.PIPE).communicate()[0]
    lines_TM = out_TM.decode()
    # print(lines_TM)

    regex_TMalign = re.compile("(?:TM-score.+)([0]\.[0-9]*)(?:.+Chain_2)")
    searchObj = re.search(regex_TMalign, lines_TM)

    # Remove useless files:
    for ext in (".sup_all_atm_lig", ".sup_all", ".sup"):
        os.remove("results/" + PU_name + ext)

    return float(searchObj.group(1))


def get_bestAlgnd_PU(nb_PU, already_selcted, peeled_pdb_name, ref_pdb_name, level):
    """
    Align the different PU (that need to be aligned) against the reference pdb
    (using TMalign) and get the number of the PU that has the maximum TMscore

    Args:
        nb_PU: Total number of PU at this given level
        already_selcted: List of the PU that have already been aligned and
        chose as best aligned PU
        peeled_pdb_name: Name (str) of the PDB that have been peeled
        ref_pdb_name: Name (str) of the PDB to align against
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
            PU_name = PEELED_PDB_ID + "_PU_" + str(level) + '_' + str(i+1)
            arr_scores[i] = TM_align(PU_name, REF_PDB_ID)

    # print(arr_scores)
    return np.argmax(arr_scores) + 1


def process_TMalign_files(bestAlgnd_PU, level, peeled_pdb_id, nb_PU):
    """
    Take the number of the best aligned PU, open the associated TM_file and with
    it:
        * Append the file gethering all aligned PU with the coordinates of this
        new PU
        * Get the resID of residues that need to be erased from the pdb (the one
        that PU are aligned against)

    Args:
        bestAlgnd_PU: The number (int) of the PU that has been determined as the
        one best aligned
        level: The current level considered (int)
        peeled_pdb_id: Name (str) of the PDB that have been peeled
        nb_PU: Total number of PU at this given level

    Returns:
        The ensemble (set) of residues that needs to be erased from the
        reference pdb (by their resID)
    """
    algnd_filename = peeled_pdb_id + '_PUs_algnd_' + str(level) + '.pdb'
    PU_name_max = peeled_pdb_id + "_PU_" + str(level) + '_' + str(bestAlgnd_PU)

    with open('results/' + PU_name_max + '.sup_atm', 'r') as sup_max:
        set_to_discard = set()

        for line in sup_max:
            if (line[0:4] == "ATOM") or ((line[0:6] == "HETATM") and
               ( (resName == "MET") or resName == "MSE") ):
               chain_ID = line[21:22].strip()

               if chain_ID == "B":
                   # Get atoms of reference pdb that are aligned, they will be
                   # discarded (corresponding to chain B):
                   set_to_discard.add(line[22:26].strip())

    with open('results/' + PU_name_max + '.sup_all_atm', 'r') as sup_all_max, \
         open('results/' + algnd_filename, 'a') as aligned_PU:

        for line in sup_all_max:
            if (line[0:4] == "ATOM") or ((line[0:6] == "HETATM") and
               ( (resName == "MET") or resName == "MSE") ):
               chain_ID = line[21:22].strip()

               if chain_ID == "A":
                   # Send chain A into the file getering all aligned PUs:
                   aligned_PU.write(line)

        aligned_PU.write("TER\n")

    # We can now delete all useless files:
    for i in range(nb_PU):
            PU_name = peeled_pdb_id + "_PU_" + str(level) + '_' + str(i+1)

            if os.path.isfile("results/" + PU_name + ".sup_atm"):
                os.remove("results/" + PU_name + ".sup_atm")
            if os.path.isfile("results/" + PU_name + ".sup_all_atm"):
                os.remove("results/" + PU_name + ".sup_all_atm")

    return set_to_discard


def erase_algned(dict_coord, set_to_discard, pdb_name):
    """
    Write a new reference pdb file, by writing only the residues that have not
    been aligned yet

    Args:
        dict_coord: Dict containing the coordinates of the reference pdb
        set_to_discard: The ensemble (set) of residues that have been aligned
        against the PU (so that must not to written)
        pdb_name: Name (str) of the reference pdb to create
    """
    with open('results/' + pdb_name + '.pdb', 'w') as pdb_file:
        resIDs = dict_coord.keys()

        for resID in resIDs:
            resID_pdb = dict_coord[resID][0]

            if resID_pdb not in set_to_discard:
                pdb_file.write(dict_coord[resID][1])


def display_curve(levels_x, TMscores_y, pdb_ref_name):
    """
    Display the curve of the values of TMscore (between aligned PU and reference
    pdb), according to the level considered

    Args:
        levels_x: Levels considered (list)
        TMscores_y: Values of TMscore (list)
        pdb_ref_name: Name (str) of the reference pdb
    """
    fig = plt.figure()
    axis = plt.subplot(111)
    plt.title("TMscore according to the level considered")
    # A list of plot to add different scores in the same plot
    # plts = []
    # for sc_name in scores:
    #     plts.append(plt.plot(tot_cumsum_res[simil_fam][sc_name]))
    plt.plot(TMscores_y, levels_x, 'bx-')
    #plt.legend([x[0] for x in plts], scores+['random'])
    plt.ylabel('TMscore with ' + pdb_ref_name)
    plt.xlabel('Level of cutting')
    axis.xaxis.set_major_locator(tck.MaxNLocator(integer=True))
    plt.show()
    #fig.savefig("results/bench/" + simil_fam + ".pdf")
    #plt.close(fig)


def calc_slope(list_x, list_y):
    x1, x2 = list_x
    y1, y2 = list_y

    return (y2-y1)/(x2-x1)

def get_best_level_2(res_levels, list_nb_PU):
    nb_levels = len(res_levels)

    if nb_levels == 1: # If only 1 level found for the protein peeled
        return 1

    else:
        i = 0
        while i < (nb_levels-1):
            slope = calc_slope(list_nb_PU[i:i+2], res_levels[i:i+2])
            i += 1


def get_best_level(res_levels, list_nb_PU):
    nb_levels = len(res_levels)

    if nb_levels == 1: # If only 1 level found for the protein peeled
        return 0

    else:
        ratio_score_nbPU = np.array(res_levels)/np.array(list_nb_PU)
        return np.argmax(ratio_score_nbPU)


def toto(ref_pdb_path, ref_pdb_id, peeled_pdb_path, peeled_pdb_id):
    """
    """
    # We need a safe copy of the ref pdb, for further TMscore use:
    os.system("cp " + ref_pdb_path + " results/" + ref_pdb_id + "_safe.pdb")

    # Creation of dssp file (needed for peeling):
    if not os.path.isfile("data/" + peeled_pdb_id + ".dss"):
        os.system("bin/dssp32 -i " + peeled_pdb_path + " > data/" +
                  peeled_pdb_id + ".dss")

    # Peeling:
    out_peel = pl.peeling(peeled_pdb_path, peeled_pdb_id)

    # Get lines from the pdb (avoid several open):
    with open(peeled_pdb_path) as pdbFile_peeled, \
         open(ref_pdb_path) as pdbFile_ref:
        dictCoord_peeled = mio.parse_pdb(pdbFile_peeled)
        dictCoord_ref = mio.parse_pdb(pdbFile_ref)


    level = 0
    res_levels = []
    list_nb_PU = []

    for line in out_peel:
        level += 1
        print("Proceeding peeling level", level)

        dict_all_PU = pl.peeled_to_dict(line)
        pl.generate_PU_pdbs(dict_all_PU, level, dictCoord_peeled, peeled_pdb_id)

        already_selcted = []
        nb_tot_PU = len(dict_all_PU)

        # Then we loop on the number of PUs, to repeat the process
        for i in range(nb_tot_PU):
        #for i in range(1):
            bestAlgnd_PU = get_bestAlgnd_PU(nb_tot_PU,
                                            already_selcted,
                                            peeled_pdb_id, ref_pdb_id,
                                            level)
            already_selcted.append(bestAlgnd_PU)

            set_to_discard = process_TMalign_files(bestAlgnd_PU, level,
                                                   peeled_pdb_id, nb_tot_PU)

            # Rewrite a pdb file (the reference one) with already
            # -aligned atoms deleted (corresponding to chain B):
            erase_algned(dictCoord_ref, set_to_discard, ref_pdb_id)

        #os.remove("PU_" + str(level) + "_algnd")
        PU_alignd_file = peeled_pdb_id + '_PUs_algnd_' + str(level) + '.pdb'
        res_levels.append(TM_score("results/" + PU_alignd_file,
                                   "results/" + ref_pdb_id + '_safe.pdb'))
        list_nb_PU.append(nb_tot_PU)


    return (res_levels, list_nb_PU)


# MAIN:
if __name__ == "__main__":
    # Get the different arguments:
    ARGS = docopt(__doc__, version='0.1')
    TO_PEELED_PDB = ARGS["--peelPdb"]
    TO_REF_PDB = ARGS["--refPdb"]

    # Creation of the results/ directory:
    if not os.path.isdir("results"):
        os.mkdir("results")

    # Extract first chain towards the results/ folder:
    # FIRST_CHAIN_ID_PEEL = "A"
    # FIRST_CHAIN_ID_REF = "A"
    PEELED_PDB, PEELED_PDB_ID = extract_chain(TO_PEELED_PDB)
    REF_PDB, REF_PDB_ID = extract_chain(TO_REF_PDB)
    PEELED_PDB_PATH = "results/" + PEELED_PDB
    REF_PDB_PATH = "results/" + REF_PDB


    # With parMATT:
    TM_parMATT = parMATT(PEELED_PDB_PATH, REF_PDB_PATH)
    print("parMATT TMscore:", TM_parMATT)

    # Peeling with TMalign:
    res_peel, list_nb_PU = toto(REF_PDB_PATH, REF_PDB_ID,
                                  PEELED_PDB_PATH, PEELED_PDB_ID)

    # Peeling with TMalign (other sense):
    res_peel_rev, list_nb_PU_rev = toto(PEELED_PDB_PATH, PEELED_PDB_ID,
                                  REF_PDB_PATH, REF_PDB_ID)

    # We start by making a simple TMalignment between both pdb:
    TMscore_ref = TM_align(PEELED_PDB_ID, REF_PDB_ID)
    for ext in ('.pdb', '.sup_atm', '.sup_all_atm'):
        os.remove("results/" + PEELED_PDB_ID + ext)
    print("Simple TMscore:", TMscore_ref)


    idx_best_level = get_best_level(res_levels, list_nb_PU)
    idx_best_level_rev = get_best_level(res_levels_rev, list_nb_PU_rev)
    print("Best peel TMscore:", res_peel[idx_best_level])
    print("Best peel TMscore (rev):", res_peel_rev[idx_best_level_rev])

    #display_curve(res_levels, range(1, len(res_levels)+1), REF_PDB_ID+" (levels)")
    display_curve(res_levels, list_nb_PU, REF_PDB_ID+" (nb_PU)")
