#!/usr/bin/env python3

"""
Module gethering all the functions that run external programs
"""

import re
import os
import subprocess as sub
import src.manage_io as mio
import urllib.request as urlreq
import pyquery as pyq


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
        The value of the associated TMscore (or -1 if no matching residues)
    """
    if peel_longer:
        cmdLine_TM = ("bin/TMscore64 " + peeled_pdb_path + " " + ref_pdb_path)
    else:
        cmdLine_TM = ("bin/TMscore64 " + ref_pdb_path + " " + peeled_pdb_path)

    out_TM = sub.Popen(cmdLine_TM.split(), stdout=sub.PIPE).communicate()[0]
    lines_TM = out_TM.decode()

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
        peeled_pdb_path: Path (str) to pdb that has been peeled (or just the
        pdb to will be moved for alignment)
        ref_pdb_path: Path (str) to the PDB to align against
        peel_longer: Boolean telling if the peeled protein is longer (or not)
        than the reference protein

    Returns:
        The value of TMscore associated with the structures aligned with
        parMATT (normalized by the longest protein)
    """
    cmdLine_parMatt = ("bin/parMATT/bin/parMatt64 " + ref_pdb_path + " " +
                       peeled_pdb_path + " -f pdb -t 1 -o output")

    out_parMatt = sub.Popen(cmdLine_parMatt.split(),
                            stdout=sub.PIPE).communicate()[0]

    # parMATT produces a single PDB file with both (aligned) structures inside
    # So we need to exctract to extract both chains (each structure) to give
    # them as argument to the TMscore program
    PEELED_PDB, PEELED_PDB_ID = mio.extract_chain("output.pdb", "A")
    REF_PDB, REF_PDB_ID = mio.extract_chain("output.pdb", "B")

    TMscore = TM_score("results/outputA.pdb", "results/outputB.pdb",
                       peel_longer)

    # Remove useless files produced by parMATT:
    os.remove("output.pdb")
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
        The value of the associated TMscore (normalized by the longest protein)
    """
    cmdLine_TM = ("bin/TMalign64 results/" + PU_name + '.pdb' +
                  " results/" + ref_pdb_name + '.pdb' + " -o " + "results/" +
                  PU_name + '.sup')

    out_TM = sub.Popen(cmdLine_TM.split(), stdout=sub.PIPE).communicate()[0]
    lines_TM = out_TM.decode()

    if peel_longer: # If peeled prot is longer, we get "normalized by chain 2"
        regex_TMalign = re.compile("(?:TM-score.+)([0]\.[0-9]*)(?:.+Chain_2)")
    else: # Else we get TMscore "normalized by chain 1"
        regex_TMalign = re.compile("(?:TM-score.+)([0]\.[0-9]*)(?:.+Chain_1)")
    searchObj = re.search(regex_TMalign, lines_TM)

    # Remove useless files:
    for ext in (".sup_all_atm_lig", ".sup_all", ".sup"):
        os.remove("results/" + PU_name + ext)

    return float(searchObj.group(1))


def gdt_pl(PU_alignd_file, ref_pdb_path, peel_longer):
    """
    Run the gdt.pl script, to get a maximized TMscore for the alignment between
    aligned PUs and reference PDB

    Args:
        PU_alignd_file: Path to the pdb file gethering all aligned PUs
        ref_pdb_path: Path (str) to the PDB to align against
        peel_longer: Boolean telling if the peeled protein is longer (or not)
        than the reference protein

    Returns:
        The value of the maximized TMscore between both structures
    """
    PU_pdb_path = "results/" + PU_alignd_file + '.pdb'
    tmp_file = PU_alignd_file + ".pdb"

    # Creation of the input file for the gdt.pl script:
    if peel_longer:
        os.system("cat " + ref_pdb_path + " > " + tmp_file)
        os.system("echo TER >> " + tmp_file)
        os.system("cat " + PU_pdb_path + " >> " + tmp_file)
    else:
        os.system("cat " + PU_pdb_path + " > " + tmp_file)
        os.system("echo TER >> " + tmp_file)
        os.system("cat " + ref_pdb_path + " >> " + tmp_file)
    os.system("echo TER >> " + tmp_file)

    cmdLine_gdt = ("perl bin/gdt.pl " + tmp_file)

    out_gdt = sub.Popen(cmdLine_gdt.split(), stdout=sub.PIPE).communicate()[0]
    lines_gdt = out_gdt.decode()

    # We get the TMscore by the chain 2, because it corresponds to the longer
    # protein
    regex_gdt = re.compile("(?:TM-score.+)([0]\.[0-9]*)(?:.+Chain 2)")
    searchObj = re.search(regex_gdt, lines_gdt)
    os.remove(tmp_file)

    return float(searchObj.group(1))


class AppURLopener(urlreq.FancyURLopener):
    version = "Mozilla/5.0"


def get_url_dom(sid_dom):
    """
    Get the exact url of the SCOP (domain) PDB, by looking through the
    HTML page associated to the domain
    """
    opener = AppURLopener()
    url_to_open = "http://scop.berkeley.edu/sid=" + sid_dom

    with opener.open(url_to_open) as response:
        pq = pyq.PyQuery(response.read())

    list_hrefs = []
    for a in pq.find('a'):
        current_href = a.get('href')
        if 'pdbstyle' in current_href and sid_dom in current_href:
            list_hrefs.append(current_href)

    if len(list_hrefs) != 1:
        print("ERROR: Several urls possible !! ")
        sys.error(2)
    else:
        return list_hrefs[0]


def dl_pdb(url_dom, pdb_id, dom_sid):
    """
    Download the given SCOP (domain) PDB file from the webpage
    """
    good_url = re.sub(r'(output=html)', 'output=txt', url_dom)

    print("Dowloading the good domain of " + pdb_id + ".pdb from the SCOP " +
          "website...")
    urlreq.urlretrieve(good_url, "data/" + dom_sid + '.pdb')
    print("Download finished !\n")
