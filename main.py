#!/usr/bin/env python3

"""PeelAlign script
Usage:
  main.py -p <peelPdb> -r <refPdb> [-b <benchMode>]

Options:
  -h --help                  help
  --version                  version of the script
  -p --peelPdb = peeled_pdb  input pdb file, that will be peeled
  -r --refPdb = ref_pdb      other input pdb file, that will be used as reference (not peeled)
  -b --benchMode = benching  Mode benchmarking (bool) [default: False]
"""


import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from docopt import docopt
import src.manage_io as mio
import src.peeling as peel
import src.external as ext


def display_curve(levels_x, levels_x_rev, TMscores_y, TMscores_y_rev,
                  ref_pdb_id, peeled_pdb_id, parMATT_TMscore, ref_TMscore):
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
    plots = []
    plots.append(plt.plot(levels_x, TMscores_y, 'bx-'))
    plots.append(plt.plot(levels_x_rev, TMscores_y_rev, 'gx-'))

    # Get maximum number of PUs:
    max_nb_PU = max(max(levels_x), max(levels_x_rev))
    # Display horizontal line for parMATT TMscore:
    plots.append(plt.plot([0, max_nb_PU],
                          [parMATT_TMscore] * 2, 'k-'))
    # Display horizontal line for reference (normal) TMscore:
    plots.append(plt.plot([0, max_nb_PU],
                          [ref_TMscore] * 2, 'm-'))

    # Add legends and axis labels:
    plt.legend([plot[0] for plot in plots],
               (peeled_pdb_id + " peeled", ref_pdb_id + " peeled",
                "parMTT TMscore", 'TMscore ref'))
    plt.ylabel('TMscores between ' + ref_pdb_id + ' and ' + peeled_pdb_id)
    plt.xlabel('Number of PUs at each level of cutting')
    axis.xaxis.set_major_locator(tck.MaxNLocator(integer=True))

    plt.show()
    #fig.savefig("results/bench/" + simil_fam + ".pdf")
    #plt.close(fig)


def check_bool_type(rep):
    """
    Check if the answer is boolean type.
    Args:
        rep: answer given by the user (str)
    Returns:
        The answer casted into a boolean
    """
    rep = rep.lower()
    if rep in ("true", "t", "yes", "y"):
        return True
    elif rep in ("false", "f", "no", "n"):
        return False
    else:
        print("Enter a boolean type (True, T, False, F) !")
        sys.exit(2)



# MAIN:
if __name__ == "__main__":
    # Get the different arguments:
    ARGS = docopt(__doc__, version='0.1')
    TO_PEELED_PDB = ARGS["--peelPdb"]
    TO_REF_PDB = ARGS["--refPdb"]
    BENCH_MODE = check_bool_type(ARGS["--benchMode"])

    # Creation of the results/ directory:
    if not os.path.isdir("results"):
        os.mkdir("results")

    # Extract first chain towards the results/ folder:
    PEELED_PDB, PEELED_PDB_ID = mio.extract_chain(TO_PEELED_PDB)
    REF_PDB, REF_PDB_ID = mio.extract_chain(TO_REF_PDB)
    PEELED_PDB_PATH = "results/" + PEELED_PDB
    REF_PDB_PATH = "results/" + REF_PDB

    # Get lines from the pdb (avoid several open):
    with open(PEELED_PDB_PATH) as pdbFile_peeled, \
         open(REF_PDB_PATH) as pdbFile_ref:
        SIZE_PEELED, DICT_COORD_PEELED = mio.parse_pdb(pdbFile_peeled)
        SIZE_REF, DICT_COORD_REF = mio.parse_pdb(pdbFile_ref)

    # Get which protein is longer than the other:
    PEEL_LONGER = False
    if SIZE_PEELED > SIZE_REF:
        PEEL_LONGER = True
    print("SIZ_PEEL:", SIZE_PEELED)
    print("SIZ_REF:", SIZE_REF)

    # With parMATT:
    TM_parMATT = ext.parMATT(PEELED_PDB_PATH, REF_PDB_PATH, PEEL_LONGER)

    # Peeled-TMalignment:
    res_peel, list_nb_PU = peel.peeled_TMalign(REF_PDB_PATH, REF_PDB_ID,
                                               DICT_COORD_REF,
                                               PEELED_PDB_PATH, PEELED_PDB_ID,
                                               DICT_COORD_PEELED, PEEL_LONGER)
    idx_best_level = peel.get_best_level(res_peel, list_nb_PU)

    # Peeling with TMalign (other sense):
    print("\nNOW REVERSE ORDER")
    tuple_res_rev = peel.peeled_TMalign(PEELED_PDB_PATH, PEELED_PDB_ID,
                                        DICT_COORD_PEELED,
                                        REF_PDB_PATH, REF_PDB_ID,
                                        DICT_COORD_REF, not PEEL_LONGER)
    res_peel_rev, list_nb_PU_rev = tuple_res_rev
    idx_best_level_rev = peel.get_best_level(res_peel_rev, list_nb_PU_rev)

    # Simple TMalignment between both pdb:
    TMscore_ref = ext.TM_align(PEELED_PDB_ID, REF_PDB_ID, PEEL_LONGER)
    # Cleaning
    for extension in ('.pdb', '.sup_atm', '.sup_all_atm'):
        os.remove("results/" + PEELED_PDB_ID + extension)
    os.remove("results/" + REF_PDB_ID + '.pdb')

    # Display the different values of scores:
    print("Simple TMscore:", TMscore_ref)
    print("Best peeled TMscore:", res_peel[idx_best_level])
    print("Best peeled TMscore (rev):", res_peel_rev[idx_best_level_rev])
    print("Max of peeled TMscores", max(res_peel[idx_best_level],
                                        res_peel_rev[idx_best_level_rev]))
    print("parMATT TMscore:", TM_parMATT, '\n')

    # Plot of the curves associated with the peeled-TMalign:
    if not BENCH_MODE:
        display_curve(list_nb_PU, list_nb_PU_rev, res_peel, res_peel_rev,
                      REF_PDB_ID, PEELED_PDB_ID, TM_parMATT, TMscore_ref)
    else:
        pass
        # Write file gethering all info for bench
