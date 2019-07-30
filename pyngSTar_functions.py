import os
import pickle
import subprocess
import pandas as pd
from Bio import SeqIO
from pyfaidx import Fasta

## Classes ##
class Allele:
	def __init__(self, gene, allele):
		self.gene = gene
		self.allele = allele

## Functions ##
def readAlleles(path):
	allelesDB = {}
	for files in os.listdir(path):
		if files.endswith('_alleles.fasta'):
			gene = files.replace('_alleles.fasta', '')
			pathfile = path+'/'+files
			with open(pathfile, 'r') as fasta:
				for record in SeqIO.parse(fasta, 'fasta'):
					name = record.name
					name = name.split('_')[1]
					if 'penA' not in record.name:
						name = name.split('.')[0]
					seq = record.seq
					rc = revComp(seq)
					allelesDB[str(seq)] = Allele(gene, name)
					allelesDB[str(rc)] = Allele(gene, name)
	return allelesDB

def readProfiles(path):
	profilesDB = {}
	pathfile = path+'/'+'ngstar_profiles.tab'
	with open(pathfile, 'r') as tab:
		for line in tab:
			linesplit = line.rstrip().split('\t')
			st = linesplit[0]
			profile = '\t'.join(linesplit[1:8])
			profilesDB[profile] = st
	return profilesDB

def revComp(seq):
	basesDic = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
	seqlist = list(seq)
	revseq = list(reversed(seqlist))
	comp = [basesDic[x] for x in revseq]
	return ''.join(comp)

def assignAllele(seq, order, allelesDB): 
	results = {}
	for i in order:
		results[i] = '-'
	for allele in allelesDB:
		if allele in seq:
			results[allelesDB[allele].gene] = allelesDB[allele].allele
	return results

def reportProfile(order, results):
	profile = results[order[0]]
	for i in order[1:]:
		profile += '\t'+results[i]
	return profile

def blastNewAlleles(query, subject, path):
	# run makeblastdb on subject (genome)
	makeblastdb_cmd = ["makeblastdb", "-in", subject, "-dbtype", "nucl", "-out", "tmp/"+subject]
	subprocess.call(makeblastdb_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	# run blastn
	blast_cmd = ["blastn", "-out", 'tmp/'+subject+'.blastn', "-outfmt", "6", "-query", path+'/'+query, "-db", "tmp/"+subject, "-evalue", "0.001"]
	subprocess.call(blast_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	# read tabular output with pandas
	bresult = pd.read_csv('tmp/'+subject+'.blastn', sep="\t", header=None)
	# get value of max bitscore
	max_bitscore = bresult[11].max()
	# get alleles with same bitscore
	subresult = bresult.loc[bresult[11] == max_bitscore]
	sublist = subresult[[0]][0].tolist()
	# get coordinates of queried alleles
	coords = subresult.iloc[0,8:10].tolist()
	# get contig name of hit
	contigloc = subresult.iloc[0,1]
	# get number from closest allele names
	clean_alleles = []
	for x in sublist:
		newx = x.split('_')[1]
		if 'penA' not in x:
			newx = newx.split('.')[0]
		clean_alleles.append(newx)
	return ['|'.join(clean_alleles[0:3])+'*', coords, contigloc]

def printNewAlleleSeqs(gene, coords, contigloc, fasta, allout, path):
	strand = 1
	startcoord, endcoord = coords
	if startcoord > endcoord:
		strand = -1
		endcoord, startcoord = coords
	# read genome and extract locus
	genome = Fasta(fasta)
	locus = genome[contigloc][startcoord:endcoord]
	if strand == -1:
		locus.seq = revComp(locus.seq)
	# print to file
	outfile = fasta.split('/').pop()+'.'+gene+'.fasta'
	if allout:
		with open(path+'/'+outfile, 'w') as out:
			out.write('>'+gene+'_'+locus.name+'_'+str(startcoord)+':'+str(endcoord)+'\n')
			out.write(locus.seq+'\n')
	return locus
