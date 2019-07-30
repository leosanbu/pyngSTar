import ahocorasick

def AC_fast(seq, order, allelesDB, allelesAC):
	results = {}
	for i in order:
		results[i] = '-'
	for end_index, (insert_order, original_value) in allelesAC.iter(seq):
		results[allelesDB[original_value].gene] = allelesDB[original_value].allele
	return results

# code to create pickle file
#allelesAC = ahocorasick.Automaton()
#for idx,key in enumerate(allelesDB):
#	allelesAC.add_word(key, (idx, key))
#allelesAC.make_automaton()
#pickle.dump((allelesDB, allelesAC), open("ngSTar_alleles_AC.pkl", "wb"))
