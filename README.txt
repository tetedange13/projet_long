===============================
   Short Python Projects

Topic: Flexible strcutural
alignment using protein peeling

AUTHOR: Felix Vandermeeren
YEAR: (dec) 2018-2019
===============================


================
EMERGENCY README
================


Emergency procedure
-------------------
Normally, you should have the fancy documentation, generated with Doxygen,
by typing in a terminal from the root directory of the project:

            xdg-open doc/html/index.html
            OR firefox doc/html/index.html


If this command does not work, and that you have Doxygen installed on your
machine, you can regenerate the documentation with:

            doxygen doc/doxy.conf

Followed by one of the above commands


If it is still not working, the essential of the documentation is below.
Just imagine that it is nice and beautiful:

Emergency description
---------------------
The project consists in implementing a program integrating proteic units into
the process of structural alignment.


Requirements
------------
python          3.6.0 or higher
numpy           1.15.2 or higher
docopt          0.6.2 or higher
pandas          0.22.0 or higher
matplotlib      3.0.2 or higher
multoprocessing 0.70.6.1 or higher



Emergency content
-----------------

PeelAlign script
  Usage:
    main.py -p <peelPdb> -r <refPdb> [-b <benchMode>] [-c <peelChain>] [-s <refChain>]

  Options:
    -h --help                  help
    --version                  version of the script
    -p --peelPdb = peeled_pdb  input pdb file, that will be peeled
    -r --refPdb = ref_pdb      other input pdb file, that will be used as reference (not peeled)
    -b --benchMode = benching  Mode benchmarking (bool) [default: False]
    -c --peelChain = pl_chain  Chain ID of the peeled PDB [default: first]
    -s --refChain = ref_chain  Chain ID of the reference PDB [default: first]

Example
-------
To run the script, with all the parameters:
    ./main.py -p data/1aoh.pdb -r data/1jlx.pdb

Some PDB files that are usable for testing can be found in the data/ folder
