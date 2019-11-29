import pickle
import ahocorasick

def updatePKL(path, allelesDB):
	allelesAC = ahocorasick.Automaton()
	for idx,key in enumerate(allelesDB):
		allelesAC.add_word(key, (idx, key))
	allelesAC.make_automaton()
	pickle.dump((allelesDB, allelesAC), open(path+'/'+'pyngSTar_alleles_AC.pkl', 'wb'))

def AC_fast(seq, order, allelesDB, allelesAC):
	results = {}
	for i in order:
		results[i] = '-'
	for end_index, (insert_order, original_value) in allelesAC.iter(seq):
		if type(results[allelesDB[original_value].gene]) is not list:
			results[allelesDB[original_value].gene] = [allelesDB[original_value].allele]
		else:
			results[allelesDB[original_value].gene].append(allelesDB[original_value].allele)
	return results

