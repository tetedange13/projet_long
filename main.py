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

import sys, os
import src.manage_io as mio
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import src.peeling as peel
import numpy as np
import src.external as ext


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
    PEELED_PDB, PEELED_PDB_ID = mio.extract_chain(TO_PEELED_PDB)
    REF_PDB, REF_PDB_ID = mio.extract_chain(TO_REF_PDB)
    PEELED_PDB_PATH = "results/" + PEELED_PDB
    REF_PDB_PATH = "results/" + REF_PDB


    # With parMATT:
    TM_parMATT = ext.parMATT(PEELED_PDB_PATH, REF_PDB_PATH)
    print("parMATT TMscore:", TM_parMATT)

    # Peeling with TMalign:
    res_peel, list_nb_PU = peel.toto(REF_PDB_PATH, REF_PDB_ID,
                                  PEELED_PDB_PATH, PEELED_PDB_ID)

    # Peeling with TMalign (other sense):
    res_peel_rev, list_nb_PU_rev = peel.toto(PEELED_PDB_PATH, PEELED_PDB_ID,
                                        REF_PDB_PATH, REF_PDB_ID)

    # We start by making a simple TMalignment between both pdb:
    TMscore_ref = ext.TM_align(PEELED_PDB_ID, REF_PDB_ID)
    for extension in ('.pdb', '.sup_atm', '.sup_all_atm'):
        os.remove("results/" + PEELED_PDB_ID + extension)
    print("Simple TMscore:", TMscore_ref)


    idx_best_level = get_best_level(res_levels, list_nb_PU)
    idx_best_level_rev = get_best_level(res_levels_rev, list_nb_PU_rev)
    print("Best peel TMscore:", res_peel[idx_best_level])
    print("Best peel TMscore (rev):", res_peel_rev[idx_best_level_rev])

    #display_curve(res_levels, range(1, len(res_levels)+1), REF_PDB_ID+" (levels)")
    display_curve(res_levels, list_nb_PU, REF_PDB_ID+" (nb_PU)")
