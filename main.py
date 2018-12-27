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
    # To view superimposed C-alpha traces of aligned regions: TM.sup
    #    To view superimposed C-alpha traces of all regions: TM.sup_all
    #    To view superimposed full-atom structures of aligned regions: TM.sup_atm
    #    To view superimposed full-atom structures of all regions: TM.sup_all_atm
    #    To view superimposed full-atom structures of all regions with ligands:
    #        TM.sup_all_atm_lig

    cmdLine_TM = ("bin/TMalign32 results/" + PU_name + '.pdb' + " results/" +
                  ref_pdb_name + '.pdb' + " -o " + "results/" + PU_name + '.sup')
    out_TM = sub.Popen(cmdLine_TM.split(), stdout=sub.PIPE).communicate()[0]
    lines_TM = out_TM.decode()
    # print(lines_TM)

    regex_TMscore = re.compile("(?:TM-score.+)([0]\.[0-9]*)(?:.+Chain_2)")
    searchObj = re.search(regex_TMscore, lines_TM)

    # Remove useless files:
    os.remove("results/" + PU_name + ".sup_all_atm_lig")
    os.remove("results/" + PU_name + ".sup_all")
    os.remove("results/" + PU_name + ".sup")

    return float(searchObj.group(1))


def TM_to_pdb(pdb_name):
    """
    Open an output file from TMalign and extract the chain A, corresponding to
    the 1st pdb now aligned on the reference pdb

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

    # Remove ".sup_all_atm" file (now useless):
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


def generate_PU_pdbs(dict_PU, level_cut, dict_coord_peeled, peeled_pdb_name):
    """
    Generate different pdb file, associated to each PU, based to the boundaries
    given as output of the peeling program

    Args:
        dict_PU: Dict containing the boundaries of each PU at a given level
        level_cut: The current level considered (int)
        dict_coord_peeled: Coordinates of the peeled pdb (dict)
        peeled_pdb_name: Name (str) of the PDB that have been peeled
    """
    nb_PU = len(dict_PU)

    for i in range(1, nb_PU+1):
        out_file = "results/" + peeled_pdb_name + "_PU_" + str(level_cut)
        with open(out_file + '_' + str(i) + '.pdb', 'w') as out_PU:
            inf_bound, sup_bound = dict_PU[i]

            for resID in range(inf_bound, sup_bound+1):
                out_PU.write(dict_coord_peeled[resID][1])


def get_bestAlgnd_PU(dict_PU, already_selcted, peeled_pdb_name, ref_pdb_name, level):
    """
    Align the different PU (that need to be aligned) against the reference pdb
    (using TMalign) and get the number of the PU that has the maximum TMscore

    Args:
        dict_PU: Dict containing the boundaries of each PU at a given level
        already_selcted: List of the PU that have already been aligned and
        chose as best aligned PU
        peeled_pdb_name: Name (str) of the PDB that have been peeled
        ref_pdb_name: Name (str) of the PDB to align against
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
            PU_name = PEELED_PDB_NAME + "_PU_" + str(level) + '_' + str(i+1)
            arr_scores[i] = TM_align(PU_name, REF_PDB_NAME)

    nb_maxScore_PU = np.argmax(arr_scores) + 1

    # We can now delete all useless files:
    for i in range(nb_PU):
        if (i+1) != nb_maxScore_PU:
            PU_name = PEELED_PDB_NAME + "_PU_" + str(level) + '_' + str(i+1)

            if os.path.isfile("results/" + PU_name + ".sup_atm"):
                os.remove("results/" + PU_name + ".sup_atm")
            if os.path.isfile("results/" + PU_name + ".sup_all_atm"):
                os.remove("results/" + PU_name + ".sup_all_atm")

    return nb_maxScore_PU


def process_TMalign_files(bestAlgnd_PU, level, peeled_pdb_name):
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
        peeled_pdb_name: Name (str) of the PDB that have been peeled

    Returns:
        The ensemble (set) of residues that needs to be erased from the
        reference pdb (by their resID)
    """
    algnd_filename = 'PU_' + str(level) + '_algnd.pdb'
    PU_name_max = peeled_pdb_name + "_PU_" + str(level) + '_' + str(bestAlgnd_PU)

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

        # aligned_PU.write("TER\n")
    os.remove('results/' + PU_name_max + '.sup_atm')
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


def TM_score(level, ref_pdb_name):
    """
    Using the TMscore binary, calculate the TMscore of a given pdb, containing
    all aligned PUs and the reference pdb (so no alignment is processed)

    Args:
        level: The current level considered (int)
        ref_pdb_name: Name (str) of the PDB to align against

    Returns:
        The value of the associated TMscore
    """
    PU_alignd_file = "PU_" + str(level) + "_algnd"  + '.pdb'
    cmdLine_TM = ("bin/TMscore32 results/" + PU_alignd_file + " results/" +
                  ref_pdb_name + '.pdb')
    out_TM = sub.Popen(cmdLine_TM.split(), stdout=sub.PIPE).communicate()[0]
    lines_TM = out_TM.decode()
    # print(lines_TM)

    regex_TMscore = re.compile("(?:TM-score.+= )([0]\.[0-9]*)")
    searchObj = re.search(regex_TMscore, lines_TM)

    return float(searchObj.group(1))


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



# MAIN:
if __name__ == "__main__":
    # Get the different arguments:
    ARGS = docopt(__doc__, version='0.1')
    PEELED_PDB_PATH = ARGS["--peelPdb"]
    REF_PDB_PATH = ARGS["--refPdb"]

    PEELED_PDB = os.path.basename(PEELED_PDB_PATH)
    REF_PDB = os.path.basename(REF_PDB_PATH)
    PEELED_PDB_NAME = os.path.splitext(PEELED_PDB)[0] #"1aoh"
    REF_PDB_NAME = os.path.splitext(REF_PDB)[0] #"1jlx"

    # Creation of dssp file:
    if not os.path.isfile("data/" + PEELED_PDB_NAME + ".dss"):
        os.system("bin/dssp " + PEELED_PDB_PATH + " > data/" + PEELED_PDB_NAME +
                  ".dss")
    if not os.path.isfile("data/" + REF_PDB_NAME + ".dss"):
        os.system("bin/dssp " + REF_PDB_PATH + " > data/" + REF_PDB_NAME +
                  ".dss")
    sys.exit()

    with open(PEELED_PDB_PATH) as pdbFile_peeled, \
         open(REF_PDB_PATH) as pdbFile_ref:
        dictCoord_1 = mio.parse_pdb(pdbFile_peeled)
        dictCoord_2 = mio.parse_pdb(pdbFile_ref)


    a_la_fac = False
    if a_la_fac:
        # Peeling:idxMax+1
        cmdLine_peel = ("bin/peeling11_4.1 -pdb " + PEELED_PDB_PATH +
                        " -dssp data/" + REF_PDB_NAME + ".dss"
                        " -R2 98 -ss2 8 -lspu 20 -mspu 0 -d0 6.0 -delta 1.5"
                        " -oss 0 -p 0 -cp 0 -npu  16")
        # The split function of the shlex module is used to generate a list of args:
        # os.system(cmd_line)
        #out, err = sub.Popen(shx.split(cmd_line), stdout=sub.PIPE).communicate()
        outPeel_1 = sub.Popen(cmdLine_peel.split(), stdout=sub.PIPE).communicate()[0]
        lines_peel = outPeel_1.decode().split('\n')

    else:
        with open("data/" + PEELED_PDB_NAME + '_peeled.txt', 'r') as outPeel_1:
            content_peel = outPeel_1.read()
            print(content_peel)
            lines_peel = content_peel.split('\n')

    # Creation of the results/ directory:
    if not os.path.isdir("results"):
        os.mkdir(results)

    # Copy pdb files towards results/ :
    os.system("cp " + PEELED_PDB_PATH + " results/")
    os.system("cp " + REF_PDB_PATH + " results/")

    # We start by making a simple TMalignment between both pdb:
    TMscore_ref = TM_align(PEELED_PDB_NAME, REF_PDB_NAME)
    TM_to_pdb(PEELED_PDB_NAME)
    print("REF", TMscore_ref)


    level = 0
    res_levels = []

    for line in lines_peel:
        if line and line[0] != '#': # There is 1 element with empty str
            level += 1

            #if level == 3:
            if True:
                dict_all_PU = peeled_to_dict(line)
                generate_PU_pdbs(dict_all_PU, level, dictCoord_1, PEELED_PDB_NAME)

                already_selcted = []
                nb_tot_PU = len(dict_all_PU)

                # Then we loop on the number of PUs, to repeat the process
                for i in range(nb_tot_PU):
                #for i in range(2):
                    bestAlgnd_PU = get_bestAlgnd_PU(dict_all_PU, already_selcted,
                                                    PEELED_PDB_NAME, REF_PDB_NAME, level)
                    already_selcted.append(bestAlgnd_PU)

                    set_to_discard = process_TMalign_files(bestAlgnd_PU, level, PEELED_PDB_NAME)

                    # Rewrite a pdb file (the reference one) with already aligned atoms
                    # deleted (corresponding to chain B):
                    erase_algned(dictCoord_2, set_to_discard, REF_PDB_NAME)

                #os.remove("PU_" + str(level) + "_algnd")
                res_levels.append(TM_score(level, REF_PDB_NAME))

    print(res_levels)
    display_curve(res_levels, range(1, len(res_levels)+1), REF_PDB_NAME)
