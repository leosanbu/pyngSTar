# pyngSTar

NG-STAR (*Neisseria gonorrhoeae* Sequence Typing for Antimicrobial Resistance) is a typing scheme published by Walter Demczuk and colleagues (J Clin Microbiol 2017 55:1454–1468. https://doi.org/10.1128/JCM.00100-17) that targets genes associated to cephalosporin, macrolides and fluoroquinolones resistance:

*penA*, *mtrR*, *porB*, *ponA*, *gyrA*, *parC* and *23S rRNA*

The NG-STAR database and web application is hosted by the Public Health Agency of Canada, National Microbiology Laboratory (walter.demczuk@canada.ca) and can be accessed from: https://ngstar.canada.ca. More information on the scheme can be found in that same webpage and their paper cited above.

**pyngSTar** is a python script to do NG-STAR typing in *Neisseria gonorrhoeae* genome assemblies. It can be run on a slow(er) or a fast mode (-f) if you have a large collection. 
* **Fast** mode (-f, recommended): uses a python module implementing the Aho-Corasick algorithm for multi-pattern string search (https://github.com/WojciechMula/pyahocorasick/).
* **Slow**(er) mode: loops through the keys of a dictionary containing the forward and reverse complementary of the alleles of the seven genes.

The main script is called **pyngSTar.py** and is accompanied by two modules:
* **pyngSTar_functions.py**: contains basic functions for reading input files, assigning allele numbers (slow mode), reporting profiles, blasting and printing new allele sequences.
* **pyngSTar_AhoCorasick.py**: imports the pyahocorasick module for fast multi-pattern string search. This is provided separately in case anyone has problems with installing the library (unlikely) and decide to use the slow mode.

### Usage:

```
usage: pyngSTar.py [options]

Get NG-STAR types from N. gonorrhoeae assemblies

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Assembly files
  -p PATH, --path PATH  Path to database containing alleles and profiles
  -f, --fast            Uses the Aho-Corasick algorithm for fast searching
                        (optional)
  -a, --alleles_out     Print fasta files with new alleles (optional)
  -o OUT_FILENAME, --out_filename OUT_FILENAME
                        Name of file to print output to (optional, default:
                        screen output)
```
| Option | Description |
| :---: | :--- |
| -i | list of assembly files (* can be used to input multiple files in a directory). |
| -p | path to a folder containing the database files. This folder must contain 9 files. |
| -f | fast allele searching option using the Aho-Corasick algorithm (optional but recommended). |
| -a | print new allele sequences if found (optional). |
| -o | output filename (default is output to screen). |

A copy of the database folder is available in pyngSTarDB.tar.gz and contains:
* \<gene\>_alleles.fasta: 7 allele sequence files with this structure. 
* ngstar_profiles.tab: 1 file containing the profiles.
* pyngSTar_alleles_AC.pkl: 1 pickle file containing the preloaded dictionary of allele sequences and numbers as well as the automaton object required by the fast searching algorithm.


### Summary of python dependencies:
* os, subprocess, argparse
* pickle
* pandas
* SeqIO
* pyfaidx
* pyahocorasick (optional)
