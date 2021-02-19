#!/usr/bin/python3
import sys
import os
import numpy as np
import csv
import pandas as pd

def get_pdbs(actual_directory,avp_directory, feature_directory):
	for file in os.listdir(actual_directory):
		filename=os.path.basename(file)
		print(filename)
		fileo=open(os.path.join(actual_directory,file))
		if filename.endswith(".pdb")== False:
			print("this file is not in correct format (i.e. not PDB){}".format(filename))
			continue
		copy=False
		seq=[]
		for line in fileo:
			if "ATOM" in line:
				if "H  95" in line:
					copy = True
					seq.append(line)
					continue
				elif "H 102" in line:
					seq.append(line)
					copy = False
					continue
				elif copy:
					seq.append(line)
		#print(seq)
		file_loop = open(os.path.join(feature_directory,"voids","file_loop",filename+'.pdb'), "a")
		file_loop.writelines(seq)
		file_loop.close()
		#bla=input("first")
		#get_voids(feature_directory,filename,file,actual_directory)

	
def get_voids(feature_directory):
	#cwd=os.getcwd()
	#print(cwd)
	#os.chdir(avp_directory)
	for file in os.listdir(os.path.join(feature_directory,"voids","file_loop")):
		filename=os.path.basename(file)
		command=os.popen("avp -R -g 0.5 -p 1.4 -n {} {} {}".format(os.path.join(feature_directory,"voids","file_loop","atom_loop",filename+".atoms_loop"),os.path.join(feature_directory,"voids","file_loop",file),os.path.join(feature_directory,"voids","file_loop","voids_loop",filename+".voids"))).readlines()
	#fileo=open(os.path.join(feature_directory,"file_loop"+'.pdb'))
	#bla=input("HEY")




avp_directory=str("/home/lilian/Desktop/avp/avp_V1.5/src/")
get_pdbs(sys.argv[1],avp_directory, sys.argv[2])
get_voids(sys.argv[2])