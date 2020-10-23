#!/usr/bin/env python
#do chmod +x name_of_script before running


import os
import sys

light=[]
heavy=[]

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
	global string
	for i in zip(position, residue, chain_type):
		one_res=one_letter_code(residue)
		string+= '{}{}\t{}\n'.format(chain_type,position,one_res)
	return (string)

for filename in os.listdir(sys.argv[1]):
	string=[]
	file=open(sys.argv[1]+"/"+filename)
	if filename.endswith(".pdb")== False:
		print("this file is not in correct format (i.e. not PDB){}".format(filename))
		continue
	for line in file:
		if "ATOM" in line and "CA" in line:
			fields = line.strip().split()
			pos= fields[5]
			res= fields[3]
			chain= fields [4]
			string= parser(pos, res, chain)
	file.close()

	new_input_file= open("{}/{}.seq".format(sys.argv[2], filename[:-4]),"w")
	for i in string:
		new_input_file.write(i)
