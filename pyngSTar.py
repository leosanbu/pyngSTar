import os
import sys
import pickle
import subprocess
import argparse as arg
from pyngSTar_functions import *
from pyngSTar_AhoCorasick import *

parser = arg.ArgumentParser(description='Get NG-STAR types from N. gonorrhoeae assemblies', usage = '%(prog)s [options]')
parser.add_argument('-i', '--input', help='Assembly files', required=False, nargs='+') # add multiple assemblies
parser.add_argument('-p', '--path', help='Path to database containing alleles and profiles', required=True)
parser.add_argument('-f', '--fast', help='Uses the Aho-Corasick algorithm for fast searching (optional, default: False)', required=False, action='store_true')
parser.add_argument('-a', '--alleles_out', help='Print fasta files with new alleles (optional, default: False)', required=False, action='store_true')
parser.add_argument('-o', '--out_filename', help='Name of file to print output to (optional, default: screen output)', required=False)
parser.add_argument('-u', '--update', help='Re-create .pkl file after updating the database', required=False, action='store_true')
arg = parser.parse_args()

# Get arguments #
filelist = arg.input
out_path = os.getcwd()
db_path = arg.path

if arg.input:
	filepath = '/'.join(filelist[0].split('/')[:-1])
	if filepath == '':
		filepath = out_path # files are in the working directory

ac = False
if arg.fast:
	ac = True

allout = False
if arg.alleles_out:
	allout = True

outfile = False
if arg.out_filename:
	outfile = arg.out_filename

update = False
if arg.update:
	update = True

## Main script ##

# read profiles file
profilesDB = readProfiles(db_path)

if update:
	allelesDB = readAlleles(db_path)
	updatePKL(db_path, allelesDB)
	print('Database (pickle file) updated in '+db_path+'/'+'pyngSTar_alleles_AC.pkl')
	sys.exit()
else:
	# Load pickled database - dictionary and automaton
	allelesDB, allelesAC = pickle.load(open(db_path+'/'+'pyngSTar_alleles_AC.pkl', 'rb'))

# order of loci in NG-STAR scheme
order = ['penA', 'mtrR', 'porB', 'ponA', 'gyrA', 'parC', '23S']

# create and print output header
header = ['strain', 'ngSTar']+order
if outfile:
	outfilehandle = open(out_path+'/'+outfile, 'w+')
	outfilehandle.write('\t'.join(header)+'\n')
else:
	print('\t'.join(header))

# create tmp/ folder for intermediate operations
subprocess.call(['mkdir', 'tmp'])

# process files
need_blast = 0
for f in filelist:
	fname = f.split('/').pop()
	with open(f, 'r') as fasta:
		concat = ''
		for record in SeqIO.parse(fasta, 'fasta'):
			concat += record.seq
		if ac:
			resultsDB = AC_fast(str(concat), order, allelesDB, allelesAC)
		else:
			resultsDB = assignAllele(concat, order, allelesDB)
		# run Blast to get closest alleles to missing loci
		save_new_alleles = {}
		for x in resultsDB:
			if resultsDB[x] == '-':
				query_file = x+'_alleles.fasta'
				subject_file = f
				blastout = blastNewAlleles(query_file, subject_file, db_path)
				if len(blastout)>1:
					need_blast = 1
					closest_allele, coords, contigloc = blastout
					resultsDB[x] = closest_allele
					save_new_alleles[x] = printNewAlleleSeqs(x, coords, contigloc, f, allout, out_path)
		# get profile and assign ST
		ngstar_prof = reportProfile(order, resultsDB)
		defaultST = 'NEW'
		ngstar_st = profilesDB.get(ngstar_prof, defaultST)
		if '-' in ngstar_prof:
			ngstar_st = '-'
		if outfile:
			outfilehandle.write(fname+'\t'+ngstar_st+'\t'+ngstar_prof+'\n')
		else:
			print(fname+'\t'+ngstar_st+'\t'+ngstar_prof)

# clean output files
subprocess.call(['rm', '-r', 'tmp'])
if need_blast == 1:
	subprocess.call(''.join('rm '+filepath+'/'+'*.fai'), shell=True)

