"""/*************************************************************************
   Program:    Qualiloop Library
   File:       myfunctions.py
   
   Version:    V1.5.1
   Date:       07.01.21
   Function:   General type definitions, defines, macros and globals
   
   Copyright:  (c) UCL, Prof. Andrew C. R. Martin 1994-2021
   Author:     Lilian Denzler
   Address:    Biomolecular Structure & Modelling Unit,
			   Department of Biochemistry & Molecular Biology,
			   University College,
			   Gower Street,
			   London.
			   WC1E 6BT.
   EMail:      andrew@bioinf.org.uk, lilian.denzler.17@ucl.ac.uk
			   
**************************************************************************
   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC
   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 
   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.
**************************************************************************
   Description:
   ============
   The library contains all functions needed to extract features of the input 
   PDB model file and to predict the RMSD score of the model using a machine 
   learning model. 
**************************************************************************
   Usage:
   ======
   This library is intended for the Qualiloop application. 
**************************************************************************
   Revision History:
   =================
   V0.1   21.01.21 Original
*************************************************************************/"""
import sys
import os
import pandas as pd
import numpy as np


#import save_RMS_lib
from qualiloop import get_abYmod2_lib
from qualiloop import seq_parser_lib
from qualiloop import charge_lib
from qualiloop import happiness_lib
from qualiloop import accessibility_lib

import doctest

"""************
	This function takes the user input PDB model file and creates a dataframe of the sequence data. 

	###change to include more of the PDB values
	***************"""
def get_input_seq(file):
	"""# return df of input sequence data
	>>> file=os.path.join("/home/lilian/sync_project/CDRH3lib/tests/1F11_1.pdb")
	>>> get_input_seq(str(file)).shape #doctest 
	(232, 3)
	"""
	file=str(file)
	input_seq_df=seq_parser_lib.run(file)
	return (input_seq_df)



def get_abymod_input(upload_dir, input_seq_df):
	"""# return csv of input sequence for abymod
	>>> current="/home/lilian/sync_project/CDRH3lib/tests"
	>>> upload_dir="/home/lilian/sync_project/CDRH3lib/tests"
	>>> test_file=os.path.join(current,"1F11_1.model")
	>>> input_seq_df=get_input_seq(test_file)
	>>> get_abymod_input(upload_dir, input_seq_df)#doctest: +NORMALIZE_WHITESPACE 
	'/home/lilian/sync_project/CDRH3lib/tests/abymod_input.csv'
	"""
	input_seq_df["position"] = input_seq_df["chain"] + input_seq_df["pos"]
	abymod_input_df=input_seq_df[['position', 'res']].copy()
	outfile = str(os.path.join(upload_dir, "abymod_input.csv"))
	abymod_input_df.to_csv(outfile, index = False, header = False, sep = '	', encoding = 'utf-8')
	return (str(outfile))



def get_models (abymod_input_csv,upload_dir, ID):
	"""# return csv of input sequence for abymod
	>>> current="/home/lilian/sync_project/CDRH3lib/tests"
	>>> upload_dir="/home/lilian/sync_project/CDRH3lib/tests"
	>>> ID ="test_model"
	>>> test_file=os.path.join(current,"1F11_1.model")
	>>> input_seq_df=get_input_seq(test_file)
	>>> abymod_input_path= get_abymod_input(upload_dir, input_seq_df)
	>>> get_models(abymod_input_path,upload_dir, ID)
	('/home/lilian/sync_project/CDRH3lib/tests/test_model.pdb.model', '/home/lilian/sync_project/CDRH3lib/tests/test_model.log')
	"""
	(model_path,log_path)=get_abYmod2_lib.pass_commands(abymod_input_csv, upload_dir, ID)
	return(model_path,log_path)




"""************
	Outputs df of sequence of input PDB model file. 
	***************"""


def get_loop_seq (input_seq_df):
	"""test
	>>> current="/home/lilian/sync_project/CDRH3lib/tests"
	>>> test_file=os.path.join(current,"1F11_1.model")
	>>> input_seq_df=get_input_seq(test_file)
	>>> get_loop_seq(input_seq_df) # doctest: +NORMALIZE_WHITESPACE
		pos res
	0    95   D
	1    96   Y
	2    97   G
	3    98   S
	4    99   T
	5   100   Y
	6  100A   G
	7  100B   F
	8   101   A
	9   102   Y
	"""
	columns=["pos","res"]
	copy = False
	seq=[]
	pos=[]
	for index, row in input_seq_df.iterrows():
		if row["chain"]=="H" and row["pos"]=="95":
			copy = True
			seq+=row["res"] 
			pos+=[row["pos"]]
			continue
		elif row["chain"]=="H" and row["pos"]=="102":
			seq+=row["res"] 
			pos+=[row["pos"]]
			copy = False
			continue
		elif copy:
			seq+=row["res"] 
			pos+=[row["pos"]]
	loop_seq_df=pd.DataFrame({'pos':pos, 'res':seq},columns=columns, index=None)
	seq_string = loop_seq_df['res'].tolist()
	seq_string="".join(seq_string)
	seq_string=pd.DataFrame([[seq_string]],columns=["seq"], index=None)

	return (loop_seq_df, seq_string)



"""************
	Outputs df of loop length and ID.
	***************"""
def get_loop_length (loop_seq_df):
	"""test
	>>> current="/home/lilian/sync_project/CDRH3lib/tests"
	>>> test_file=os.path.join(current,"1F11_1.model")
	>>> input_seq_df=get_input_seq(test_file)
	>>> loop_seq_df, seq_string=get_loop_seq(input_seq_df) 
	>>> get_loop_length(loop_seq_df) # doctest: +NORMALIZE_WHITESPACE
	length
	0      10
	"""
	columns=["length"]
	length=loop_seq_df.shape[0]
	length_df = pd.DataFrame([length], columns=columns)
	return (length_df)



"""************
	Outputs df of loop total charge, nr of charged residues in the loop.
	***************"""

def get_loop_charge(loop_seq_df):
	"""test
	>>> loop_seq_df= pd.DataFrame(np.array([["95","D"],["96","Y"],["97","G"],["98","S"],["99","T"],["100","Y"],["100A","G"],["100B","F"],["101","A"],["102","Y"]]), columns=['pos', 'res'])
	>>> get_loop_charge(loop_seq_df) # doctest: +NORMALIZE_WHITESPACE
	total_charge  nr_charged
	0            -1           1
	"""
	charge_df=charge_lib.get_charge(loop_seq_df)
	return (charge_df)


"""************
	Outputs df of Happiness_score
	***************"""
'''
def get_happiness_score(file):
	#only works on ACRM server 
	"""test
	>>> loop_seq_df= pd.DataFrame(np.array([["95","D"],["96","Y"],["97","G"],["98","S"],["99","T"],["100","Y"],["100A","G"],["100B","F"],["101","A"],["102","Y"]]), columns=['pos', 'res'])
	>>> get_happiness_score(loop_seq_df) # doctest: +NORMALIZE_WHITESPACE
	0
	"""
	# You must set DATADIR to the BiopLib/BiopTools data directory  ??
	######################
	happiness_df=happiness_lib.get_happy(file)
	return (happiness_df)

'''



"""************
	Outputs df of identity and similarity to template sequence. 
	***************"""
def similarity_score(log_name,upload_dir):
	"""test
	>>> current="/home/lilian/sync_project/CDRH3lib/tests"
	>>> upload_dir="/home/lilian/sync_project/CDRH3lib/tests"
	>>> ID ="test_model"
	>>> test_file=os.path.join(current,"1F11_1.model")
	>>> input_seq_df=get_input_seq(test_file)
	>>> abymod_input_path= get_abymod_input(upload_dir, input_seq_df)
	>>> (model_path, log_path)=get_models(abymod_input_path,upload_dir, ID)
	>>> similarity_score(log_path,upload_dir)
	"""
	columns=["identity","similarity", "template", "target"]
	data=[]
	with open(str(log_name), "r") as infile:
		for line in infile:
			#INFO: CDR-H3 (YEIR/YEWA) SeqID: 0.500 Similarity: 0.381
			if "INFO: CDR-H3" in line:
				output=line.split(" ")
				template_target=output[2]
				template_target=template_target.split("/")
				target=template_target[0]
				template= template_target[1]
				target=target.replace("(", "")
				template=template.replace(")", "")

				identity=output[4]
				similarity=output[6]
				similarity=similarity.replace('"', "")
				similarity=similarity.replace("\n","")
				write=[identity,similarity, template, target]
				data.append(write)

	df = pd.DataFrame(data, columns=columns)
	return (df)




"""************
	Outputs df of Hydrophobicity values based on the consensus values by Eisenberg et al. 
	(Eisenberg, et al 'Faraday Symp.Chem.Soc'17(1982)109). 
	output:
	Mean of hydrophobicity values of loop
	Sum of absolute Differences between loop and template loop
	***************"""				
def hydropathy(loop_seq_df, ID, log_name, upload_dir):
	"""test
	>>> ID="test"
	>>> current="/home/lilian/sync_project/CDRH3lib/tests"
	>>> test_file=os.path.join(current,"1F11_1.model")
	>>> log_name="/home/lilian/sync_project/CDRH3lib/tests/1F11_1.log"
	>>> upload_dir="/home/lilian/sync_project/CDRH3lib/tests"
	>>> input_seq_df=get_input_seq(test_file)
	>>> loop_seq_df, seq_string=get_loop_seq(input_seq_df)
	>>> hydropathy(loop_seq_df, ID, log_name, upload_dir)
	0
	"""
	columns=["Hydropathy","Hydropathy_diff","ID"]
	#Consensus values: Eisenberg, et al 'Faraday Symp.Chem.Soc'17(1982)109
	Hydrophathy_index = {'A': 00.250, 'R': -1.800, "N": -0.640, "D": -0.720, "C": 00.040, "Q": -0.690, "E": -0.620, "G": 00.160, "H": -0.400, "I": 00.730, "L": 00.530, "K": -1.100, "M": 00.260, "F": 00.610, "P": -0.070,
							"S": -0.260, "T": -0.180, "W": 00.370, "Y": 00.020, "V": 00.540, "X": -0.5}#-0.5 is average
	df_sim=similarity_score(log_name,upload_dir)
	name= ID
	hydro_value=0
	data=[]
	try:
		for index, row in loop_seq_df.iterrows():
			aminoacid=row['res']
			pos=row["pos"]
			hydro_value+=Hydrophathy_index[aminoacid]
			row=df_sim.iloc[[0]]
			if row.empty==True:
				hydro_diff=None
			else: 
				template=row["template"].values[0]
				target=row["target"].values[0]
				for a,b in zip(template, target):
					hydro_diff=abs(abs(Hydrophathy_index[b])-abs(Hydrophathy_index[a]))
	except:
		pass
	write=[hydro_value,hydro_diff,name]
	data.append(write)
	df = pd.DataFrame(data, columns=columns)
	return (df)


"""************
	Outputs df of accessibility 
	***************"""				   
def accessibility(file, ID, upload_dir):
	"""test
	>>> current="/home/lilian/sync_project/CDRH3lib/tests"
	>>> ID="test"
	>>> file=os.path.join(current,"1F11_1.model")
	>>> upload_dir="/home/lilian/sync_project/CDRH3lib/tests"
	>>> accessibility(file, ID, upload_dir) #doctest: +NORMALIZE_WHITESPACE
	0
	"""
	df=accessibility_lib.pdbsolv_run(file, ID, upload_dir)
	return (df)

def merge(file, ID, list_of_features, upload_dir):
	"""test
	>>> current="/home/lilian/sync_project/CDRH3lib/tests"
	>>> ID="test"
	>>> list_of_features=["loop_seq", "loop_length", "total_charge", "identity"]
	>>> file=str(os.path.join(current,"1F11_1.pdb"))
	>>> upload_dir="/home/lilian/sync_project/CDRH3lib/tests"
	>>> merge(file, ID, list_of_features, upload_dir) #doctest: +NORMALIZE_WHITESPACE
	0
	"""
	input_seq_df=get_input_seq(file)
	abymod_input_csv=get_abymod_input(upload_dir, input_seq_df)
	(model_path,log_path)=get_models (abymod_input_csv,upload_dir, ID)
	loop_seq_df, seq_string=get_loop_seq (input_seq_df)
	loop_length_df=get_loop_length (loop_seq_df)
	final_df=[]

	if "loop_seq" in list_of_features:
		final_df.append(seq_string)

	if "loop_length" in list_of_features:
		loop_length_df=get_loop_length (loop_seq_df)
		final_df.append(loop_length_df)

	if "total_charge" in list_of_features:
		loop_charge_df=get_loop_charge(loop_seq_df)
		final_df.append(loop_charge_df[["total_charge"]])

	if "nr_charged" in list_of_features:
		loop_charge_df=get_loop_charge(loop_seq_df)
		final_df.append(loop_charge_df[["nr_charged"]])

	'''if "tip_res" in list_of_features:

	if "protrusion" in list_of_features:'''

	if "template" in list_of_features:
		similarity_df=similarity_score(log_path,upload_dir)
		final_df.append(similarity_df[["template"]])

	if "identity" in list_of_features:
		similarity_df=similarity_score(log_path,upload_dir)
		final_df.append(similarity_df[["identity"]])

	if "similarity" in list_of_features:
		similarity_df=similarity_score(log_path,upload_dir)
		final_df.append(similarity_df[["similarity"]])

	if "hydrophobicity" in list_of_features:
		hydropathy_df=hydropathy(loop_seq_df, ID, log_path, upload_dir)
		final_df.append(hydropathy_df)

	if "Access" in list_of_features:
		accessibility_df=accessibility(file, ID, upload_dir)
		final_df.append(accessibility_df[["Access"]])

	if "Relaccess" in list_of_features:
		accessibility_df=accessibility(file, ID, upload_dir)
		final_df.append(accessibility_df[["Relaccess"]])

	if "Scaccess" in list_of_features:
		accessibility_df=accessibility(file, ID, upload_dir)
		final_df.append(accessibility_df[["Scaccess"]])

	if "Relscaccess" in list_of_features:
		accessibility_df=accessibility(file, ID, upload_dir)
		final_df.append(accessibility_df[["Screlacc"]])

	result = pd.concat(final_df, axis=1)
	result.to_csv(os.path.join(upload_dir,"features"+".csv"), index=False)
	result_path=os.path.join(upload_dir,"features"+".csv")
	return(result)

'''

#creation of binary and nominal features will be included in the modeling stage
def RMSD_binary(actual_directory, model_directory, feature_directory):
	columns=["local_AA_bin","local_CA_bin","global_AA_bin", "global_CA_bin", "ID"]
	data=[]
	threshold=int(input("Threshold value in ANgstrom (above-> not good model)"))
	RMSD_path=Path(input("enter FULL path of csv file containing RMSD values, press enter to calculate these"))
	if os.path.exists(RMSD_path) ==True:
		df=pd.read_csv(RMSD_path)
	else:
		input("AAAAH")
		df=get_RMSDs(actual_directory, model_directory, feature_directory)
	local_AA= df.local_AA.tolist()
	local_CA= df.local_CA.tolist()
	global_AA= df.global_AA.tolist()
	global_CA=df.global_CA.tolist()
	name=df.ID.tolist()
	for (a,b,c,d,n) in zip(local_AA, local_CA, global_AA, global_CA, name):
		a=float(a)
		b=float(b)
		c=float(c)
		d=float(d)
		if a>=threshold:
			a=1
		else:
			a=0
		if b>=threshold:
			b=1
		else:
			b=0
		if c>=threshold:
			c=1
		else:
			c=0
		if d>=threshold:
			d=1
		else:
			d=0
		write=[a,b,c,d,n]
		data.append(write)
	dfnew = pd.DataFrame(data, columns=columns)
	print("RMSD_binary:",dfnew.shape)
	return (dfnew)

def RMSD_nom(actual_directory, model_directory, feature_directory):
	columns=["local_AA_nom","local_CA_nom","global_AA_nom", "global_CA_nom", "ID"]
	data=[]
	RMSD_path=Path(input("enter FULL path of csv file containing RMSD values, press enter to calculate these"))
	if os.path.exists(RMSD_path) ==True:
		df=pd.read_csv(RMSD_path, header=0)
	else:
		input("AAAAH")
		df=get_RMSDs(actual_directory, model_directory, feature_directory)
	local_AA= df.local_AA.tolist()
	local_CA= df.local_CA.tolist()
	global_AA= df.global_AA.tolist()
	global_CA=df.global_CA.tolist()
	name=df.ID.tolist()
	threshlow=[0,2,4,6,8,10,12,14,16,18,20]
	threshup=[2,4,6,8,10,12,14,16,18,20,100]
	values=[1,2,3,4,5,6,7,8,9,10,11]
	print(len(threshlow), len(threshup), len(values))
	for (a,b,c,d,n) in zip(local_AA, local_CA, global_AA, global_CA, name):
		a=float(a)
		b=float(b)
		c=float(c)
		d=float(d)
		for (threshold_low,threshold_up,val) in list(zip(threshlow,threshup,values)):
			if type(a)==str:
				pass
			elif a>threshold_low and a<threshold_up:
				a=val
			if type(b)==str:
				pass
			elif b>threshold_low and b<threshold_up:
				b=val
			if type(c)==str:
				pass
			elif c>threshold_low and c<threshold_up:
				c=val
			if type(d)==str:
				pass
			elif d>threshold_low and d<threshold_up:
				d=val
		write=[int(a),int(b),int(c),int(d),n]
		data.append(write)
	dfnew = pd.DataFrame(data, columns=columns)
	print("RMSD_nom:",dfnew.shape)
	return (dfnew)

'''

if __name__ == '__main__':
	#doctest.testmod()
	ID="test"
	list_of_features=["loop_seq", "loop_length", "total_charge", "identity"]
	file=str(os.path.join("/home/lilian/sync_project/CDRH3lib/tests","1F11_1.pdb"))
	upload_dir="/home/lilian/sync_project/CDRH3lib/tests"
	print(merge(file, ID, list_of_features, upload_dir))
#alternatively test with:  python -m doctest -v myfunctions.py
