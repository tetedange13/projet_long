# Peeled-TMalignment


## Purpose of the program
This program aims in integrating **proteic units** (PU) into the process of structural alignment, in order to add **some flexibility**.


## Installation
- Clone our repository and open it
  ```
  $ git clone https://github.com/tetedange13/projet_long.git
  $ cd projet_long
  ```

- The required packages  are:
  ```
  python          3.6.0 or higher
  numpy           1.15.2 or higher
  docopt          0.6.2 or higher
  pandas          0.22.0 or higher
  matplotlib      3.0.2 or higher
  multoprocessing 0.70.6.1 or higher
  ```


## Usage
- The program is composed of 2 python scripts `main.py` and `bench.py`, all in the base directory.  
All modules are stored in the `src/` directory. <br>

- Show help (for the main):
```
  $ python3 main.py -h
```
```
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
```

- Run the program on the **RIPC dataset for benchmarking**
```
$ python3 bench.py
```

## Example of command for 2 pdb files:
  ```
  $ python3 main.py -p data/d1a2pa_.pdb -r data/d1afra_.pdb
  ```
This should generate a **nice plot (into the root directory), called `curve.pdf`**.


## Documentation
- To open the nice HTML documentation:
```
$ firefox doc/html/index.html
```


- To regenerate documentation using [Doxygen](https://github.com/doxygen/doxygen):
```
$ doxygen doc/doxy.conf
```

## Project author: [Félix Vandermeeren](https://github.com/tetedange13)
