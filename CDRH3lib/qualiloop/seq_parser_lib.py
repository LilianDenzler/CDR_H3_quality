#!/usr/bin/env python
#do chmod +x name_of_script before running

import pandas as pd

#argv[1]= direcoty with actual PDBs, argv[2] = directory to keep input sequences

def one_letter_code(residue):
	dic = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','XAA': 'X', 'UNK':'X'}
	#one_letter= atomium.structures.code()    #####doesnt work -> why?
	if len(residue) % 3 != 0:
		print(residue)
		raise ValueError("error")
	one_letter= dic[residue]
	return (one_letter)

def parser(position,residue,chain_type):
	for i in zip(position, residue, chain_type):
		one_res=one_letter_code(residue)
		if chain_type=="H" or chain_type=="L":
			string=(chain_type,position,one_res)
			return (chain_type,position,one_res)

def run(file):
	dic_input= {'chain': [],'pos':[],'res':[]}
	with open(file, "r") as file:
		for line in file:
			if "ATOM" in line and "CA" in line:
				fields = line.strip().split()
				pos= fields[5]
				res= fields[3]
				chain= fields [4]
				(chain_type,position,one_res)= parser(pos, res, chain)
				dic_input['chain'].append(chain_type)
				dic_input['pos'].append(position)
				dic_input['res'].append(one_res)
	df = pd.DataFrame(dic_input, columns = ['chain', 'pos','res'], index=None)
	return (df)

