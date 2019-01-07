#!/usr/bin/env python3

"""
Benchmarking script

Process to the benchmarking if needed and then generate nice plots and other
interesting information about the results of the benchmarking
"""
# Usage:
#   bench.py -p <peelPdb> -r <refPdb> [-b <benchMode>]
#
# Options:
#   -h --help                  help
#   --version                  version of the script
#   -p --peelPdb = peeled_pdb  input pdb file, that will be peeled
#   -r --refPdb = ref_pdb      other input pdb file, that will be used as reference (not peeled)
#   -b --benchMode = benching  Mode benchmarking (bool) [default: False]
# """

import sys
import os
import re
import numpy as np
import pandas as pd
import subprocess as sub
import matplotlib.pyplot as plt
import src.external as ext


def disp_barplot(all_means):
    fig = plt.figure()
    axis = plt.subplot(111)
    plt.title("Mean values of TMscore according to the method considered")


    barlist = axis.bar(all_means.index, all_means.values)
    fig.subplots_adjust(bottom=0.2)
    plt.xticks(rotation=45)
    for bar in barlist:
        bar.set_color(list(np.random.random(size=3)))

    fig.savefig("barplot.pdf")
    plt.close(fig)


# MAIN:
if __name__ == "__main__":
    if not os.path.isfile('bench.csv'):
        # Get dataset table:
        RIPC_txt = pd.read_csv('data/RIPC_dataset.txt', sep='\t')

        # Write header of the file containing the results of the bench:
        with open('bench.csv', 'w') as bench_file:
            bench_file.write("peel_pdb_id-ref_pdb_id;TMscore_ref;TM_parMATT;" +
                       "best_peel_TM_rev;best_peel_TM;max_peel_TM;best_peel_gdt;" +
                       "best_peel_gdt_rev;max_peel_gdt\n")


        for idx, row in RIPC_txt.iterrows():
            len_dom1, len_dom2 = row['Length1'],row['Length2']
            dom1_sid, dom2_sid = row['Domain1'], row['Domain2']
            pdb_id_dom1, pdb_id_dom2 = dom1_sid[1:5], dom2_sid[1:5]
            chainID_sid1, chainID_sid2 = dom1_sid[5], dom2_sid[5]

            # If the PDB is absent from the data/ folder, it is downloaded:
            if not os.path.isfile("data/" + dom1_sid + '.pdb'):
                url_dom1 = ext.get_url_dom(dom1_sid)
                ext.dl_pdb(url_dom1, pdb_id_dom1, dom1_sid)
            if not os.path.isfile("data/" + dom2_sid + '.pdb'):
                url_dom2 = ext.get_url_dom(dom2_sid)
                ext.dl_pdb(url_dom2, pdb_id_dom2, dom2_sid)

            cmd_main = ("./main.py -p data/" + dom1_sid + '.pdb -r data/' +
                        dom2_sid + ".pdb -b t")
            os.system(cmd_main)


    else: # Display plots and information about results:
        results = pd.read_csv('bench.csv', sep=';', index_col=0)
        all_means = results.loc[:, ['TMscore_ref', 'TM_parMATT', 'max_peel_TM',
                               'max_peel_gdt']].mean()
        disp_barplot(all_means)

        print("Couples dont le peeled-TMscore est meilleur que celui de reference:")
        TM_sup_ref = results.index[results['max_peel_TM'] > results['TMscore_ref']]
        print(results.loc[TM_sup_ref, ['TMscore_ref', 'max_peel_TM']])

        print('\n')
        print("Couples dont le peeled-gdt-TMscore est meilleur que celui de " +
              "reference:")
        gdt_sup_ref = results.index[results['max_peel_gdt'] > results['TMscore_ref']]
        print(results.loc[gdt_sup_ref, ['TMscore_ref', 'max_peel_gdt']])


        # Couples with high TMscores, according the different methods:
        cutoff_high = 0.5
        print('\n\n')

        print("Couples dont le peeled-TMscore est supérieur à " +
              str(cutoff_high) + ":")
        TM_high = results.index[results['max_peel_TM'] > cutoff_high]
        print(results.loc[TM_high, ['TMscore_ref', 'max_peel_gdt',
                                        'max_peel_TM']])

        print('\n')
        print("Couples dont le peeled-gdt-TMscore est supérieur à " +
              str(cutoff_high) + ":")
        gdt_high = results.index[results['max_peel_gdt'] > cutoff_high]
        print(results.loc[gdt_high, ['TMscore_ref', 'max_peel_gdt',
                                        'max_peel_TM']])

        print('\n')
        print("Couples dont le TMscore de reference est supérieur à " +
              str(cutoff_high) + ":")
        TM_ref_high = results.index[results['TMscore_ref'] > cutoff_high]
        print(results.loc[TM_ref_high, ['TMscore_ref', 'max_peel_gdt',
                                        'max_peel_TM']])

        # Couples with low TMscores, according the different methods:
        cutoff_low = 0.17
        print('\n\n')

        print("Couples dont le peeled-TMscore est inferieur à " + str(cutoff_low) +
              ":")
        TM_low = results.index[results['max_peel_TM'] < cutoff_low]
        print(results.loc[TM_low, ['TMscore_ref', 'max_peel_gdt',
                                        'max_peel_TM']])

        print("NB mauvais peeled-TMscore:", len(TM_low))

        print('\n')
        print("Couples dont le peeled-gdt-TMscore est inferieur à " +
              str(cutoff_low) + ":")
        gdt_low = results.index[results['max_peel_gdt'] < cutoff_low]
        print(results.loc[gdt_low, ['TMscore_ref', 'max_peel_gdt',
                                        'max_peel_TM']])
        print("NB mauvais peeled-gdt-TMscore:", len(gdt_low))

        print('\n')
        print("Couples dont le TMscore de reference est inferieur à " +
              str(cutoff_low) + ":")
        TM_ref_low = results.index[results['TMscore_ref'] < cutoff_low]
        print(results.loc[TM_ref_low, ['TMscore_ref', 'max_peel_gdt',
                                        'max_peel_TM']], '\n')
