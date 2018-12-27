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

            doxygen doxy.conf

Followed by one of the above commands


If it is still not working, the essential of the documentation is below.
Just imagine that it is nice and beautiful:

Emergency description
---------------------
The project consists in reimplementing an existing algorithm that finds
transmembrane parts from transmembrane proteins
From this article:



Requirements
------------
docopt
numpy
matplotlib



Emergency content
-----------------

Usage:
  main.py -p <peeledPdb> -r <refPdb>

Options:
  -h --help                  help
  --version                  version of the script
  -p --peeledPdb             input pdb file, that will be peeled
  -r --refPdb                other input pdb file, that will be used as reference (not peeled)

Example
-------
To run the script, with all the parameters:
    ./main.py -p data/1aoh.pdb -r data/1jlx.pdb
