#!/usr/bin/env python3

"""
Module gethering all the functions that run external programs
"""

import re
import os
import subprocess as sub
import src.manage_io as mio


def TM_score(peeled_pdb_path, ref_pdb_path, peel_longer):
    """
    Using the TMscore binary, calculate the TMscore of a given pdb, containing
    all aligned PUs and the reference pdb (so no alignment is processed)

    Args:
        peeled_pdb_path: Path (str) to pdb that has been peeled (or just the
        pdb to will be moved for alignment)
        ref_pdb_path: Path (str) to the PDB to align against
        peel_longer: Boolean telling if the peeled protein is longer (or not)
        than the reference protein

    Returns:
        The value of the associated TMscore
    """
    if peel_longer:
        cmdLine_TM = ("bin/TMscore32 " + peeled_pdb_path + " " + ref_pdb_path)
    else:
        cmdLine_TM = ("bin/TMscore32 " + ref_pdb_path + " " + peeled_pdb_path)
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


def parMATT(peeled_pdb_path, ref_pdb_path, peel_longer):
    """
    Args:
        peeled_pdb_path:
        ref_pdb_path:
        peel_longer: Boolean telling if the peeled protein is longer (or not)
        than the reference protein
    """
    cmdLine_parMatt = ("bin/parMATT/bin/parMatt32 " + peeled_pdb_path + " " +
                       ref_pdb_path + " -t 1 -o output")

    out_parMatt = sub.Popen(cmdLine_parMatt.split(),
                            stdout=sub.PIPE).communicate()[0]
    # print(out_parMatt.decode())

    PEELED_PDB, PEELED_PDB_ID = mio.extract_chain("output.pdb", "A")
    REF_PDB, REF_PDB_ID = mio.extract_chain("output.pdb", "B")

    # Remove useless files produced by parMATT:
    for extension in ('spt', 'fasta', 'txt'):
        os.remove("output." + extension)

    TMscore = TM_score("results/" + PEELED_PDB, "results/" + REF_PDB, peel_longer)
    os.remove("results/outputA.pdb")
    os.remove("results/outputB.pdb")

    return TMscore


def TM_align(PU_name, ref_pdb_name, peel_longer):
    """
    Using the TMalign binary, proceed to a structural alignment between a given
    PU and the reference pdb (TMscore is returned as the stdout of TMalign
    program)

    Args:
        PU_name: Name (str) of the PU to align
        ref_pdb_name: Name (str) of the PDB to align against
        peel_longer: Boolean telling if the peeled protein is longer (or not)
        than the reference protein

    Returns:
        The value of the associated TMscore
    """
    if peel_longer: # If peeled prot is longer, we keep the order as is
        cmdLine_TM = ("bin/TMalign32 results/" + PU_name + '.pdb' +
                      " results/" + ref_pdb_name + '.pdb' + " -o " + "results/" +
                      PU_name + '.sup')

    else: # Else we invert both proteins
        cmdLine_TM = ("bin/TMalign32 results/" + ref_pdb_name + '.pdb' +
                      " results/" + PU_name + '.pdb' + " -o " + "results/" +
                      PU_name + '.sup')

    out_TM = sub.Popen(cmdLine_TM.split(), stdout=sub.PIPE).communicate()[0]
    lines_TM = out_TM.decode()
    # print(lines_TM)

    regex_TMalign = re.compile("(?:TM-score.+)([0]\.[0-9]*)(?:.+Chain_1)")
    searchObj = re.search(regex_TMalign, lines_TM)

    # Remove useless files:
    for ext in (".sup_all_atm_lig", ".sup_all", ".sup"):
        os.remove("results/" + PU_name + ext)

    return float(searchObj.group(1))


def gdt_pl(PU_pdb_path, ref_pdb_path):
    """
    """
    os.system("cat " + PU_pdb_path + " > toto.pdb")
    os.system("echo TER >> toto.pdb")
    os.system("cat " + ref_pdb_path + " >> toto.pdb")
    os.system("echo TER >> toto.pdb")

    cmdLine_gdt = ("perl bin/gdt.pl toto.pdb")

    out_gdt = sub.Popen(cmdLine_gdt.split(), stdout=sub.PIPE).communicate()[0]
    lines_gdt = out_gdt.decode()
    print(lines_gdt)

    # We get the TMscore by the chain 2, because it corresponds to the ref PDB
    regex_gdt = re.compile("(?:TM-score.+)([0]\.[0-9]*)(?:.+Chain 2)")
    searchObj = re.search(regex_gdt, lines_gdt)

    return float(searchObj.group(1))
