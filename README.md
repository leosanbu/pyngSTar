# pyngSTar

NG-STAR (*Neisseria gonorrhoeae* Sequence Typing for Antimicrobial Resistance) is a typing scheme published by Walter Demczuk and colleagues (J Clin Microbiol 2017 55:1454–1468. https://doi.org/10.1128/JCM.00100-17) that targets genes associated to cephalosporin, macrolides and fluoroquinolones resistance:

*penA*, *mtrR*, *porB*, *ponA*, *gyrA*, *parC* and *23S rRNA*

The NG-STAR database and web application is hosted by the Public Health Agency of Canada, National Microbiology Laboratory (walter.demczuk@canada.ca) and can be accessed from: https://ngstar.canada.ca. More information on the scheme can be found in that same webpage and their paper cited above.

**pyngSTar** is a python3 script to do NG-STAR typing in *Neisseria gonorrhoeae* genome assemblies. To look for exact matches, it can be run on a slow(er) or a fast mode (-f) if you have a large collection. 
* **Fast** mode (-f, recommended): uses a python module implementing the Aho-Corasick algorithm for multi-pattern string search (https://github.com/WojciechMula/pyahocorasick/).
* **Slow**(er) mode: loops through the keys of a dictionary containing the forward and reverse complementary of the alleles of the seven genes.

The main script is called **pyngSTar.py** and is accompanied by two modules:
* **pyngSTar_functions.py**: contains basic functions for reading input files, assigning allele numbers (slow mode), reporting profiles, blasting and printing new allele sequences.
* **pyngSTar_AhoCorasick.py**: imports the pyahocorasick module for fast multi-pattern string search. This is provided separately in case anyone has problems with installing the library (unlikely) and decide to use the slow mode.

If an exact match is not found, **blastn** is called to report the closest matches and the new allele is printed to an output file if requested with -a. Closest matches are marked in the profiles table with an asterisk. If more than 3 alleles are equally close, only 3 are shown in the table. In this case, something may be wrong, i.e. the whole sequence of the allele is not found and only part of it is beign considered, so it is good to take a look at the output fasta files of potential new alleles.

Profiles not found in the local database are marked in the table as 'NEW'.

### Summary of python dependencies:

* os, subprocess, argparse
* pickle
* pandas
* SeqIO
* pyfaidx
* pyahocorasick (optional)

### Options:

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

A copy of the database (downloaded on July 2019 from https://ngstar.canada.ca) folder is available in pyngSTarDB.tar.gz and contains:
* \<gene\>_alleles.fasta: 7 allele sequence files with this structure. 
* ngstar_profiles.tab: 1 file containing the profiles.
* pyngSTar_alleles_AC.pkl: 1 pickle file containing the preloaded dictionary of allele sequences and numbers as well as the automaton object required by the fast searching algorithm.

### Usage:

To run pyngSTar in a fast mode (-f) and requesting new alleles to be printed as fasta files (-a):
```
python pyngSTar.py -f -a -i /path/to/assemblies/*.fa -p pyngSTarDB/
```
Add -o \<filename\> if you want the output to go to a file instead of the screen.

### Output example:

```
strain         ngSTar  penA     mtrR    porB  ponA  gyrA  parC  23S
10356_1#33.fa  NEW     19.001   25|86*  100   1     100   7     100
10625_6#40.fa  NEW     14.001*  54      100   1     1     18    100
10625_6#8.fa   127     5.002    1       8     1     1     3     100
11792_4#61.fa  356     2.001    10      100   100   1     18    100
```
This will produce two fasta files containing the new alleles found for *penA* and *mtrR*. The header of the sequences indicate have the gene name, the contig name where the gene was found and the start and end coordinates in that contig where the sequence was pulled out from:
```
>mtrR_.10356_1_33.8_29704:30403
TTGCATGGTTACAAAGTCTTTTTTATAATCCGCCCTCATCAAACCGACCCGAAACGAAACCGCCA...
````

For the 419 *N. gonorrhoeae* strains from Sánchez-Busó *et al.* 2019 (https://www.nature.com/articles/s41564-019-0501-y), pyngSTar on the fast mode took just 1min 31s!!

### Why another typing script?

Same as NG-MASTER (https://doi.org/10.1099/mgen.0.000076), each typing scheme has its own particularities!!

### Contact:

If you use it and it crashes, or you have an idea to improve it, i'd love to hear from you (leonor.sanchezbuso@bdi.ox.ac.uk).

### Citation:

Please, if you use this script for a publication, cite the original work by Walter Demczuk and colleagues https://doi.org/10.1128/JCM.00100-17 and a link to this github page! 
