#!/usr/bin/env python
import sys
import os
import pandas as pd


def pdbsolv_run(file, ID, upload_dir):
	columns=["Access","Relacc","Scacc","Screlacc","Access_avg","Relacc_avg","Scacc_avg","Screlacc_avg","ID"]
	data=[]
	filename=file
	nr_sad=0
	name= ID
	try:
		command=os.popen("pdbsolv -r stdout {}".format(os.path.join(file))).readlines()# RESIDUE  AA   ACCESS  RELACC  SCACC   SCRELACC
		copy=False
		results=False
		access=0
		relacc=0
		scacc=0
		screlacc=0
		counter=0
		for i in command:
			i=i.split()
			#print(i)
			if "END" in i:
				results=True
				continue
			if results==True:
				if "95"== i[2] and "H"==i[1]:
					copy = True
				if copy==True:
					counter+=1
					access+=float(i[4])
					relacc+=float(i[5])
					scacc+=float(i[6])
					screlacc+=float(i[7])
				if "102"==i[2] and "H"==i[1]:
					copy = False
		if counter==0:
			access_avg=None
			relacc_avg=None
			scacc_avg=None
			screlacc_avg=None
			access=None
			relacc=None
			scacc=None
			screlacc=None
			write=[access,relacc,scacc,screlacc, access_avg,relacc_avg, scacc_avg, screlacc_avg, name]
			data.append(write)

		else:
			access_avg=access/counter
			relacc_avg=relacc/counter
			scacc_avg=scacc/counter
			screlacc_avg=screlacc/counter

			write=[access,relacc,scacc,screlacc, access_avg,relacc_avg, scacc_avg, screlacc_avg, name]
			data.append(write)
	except:
		print("accessibility:error", name)
			
	df = pd.DataFrame(data, columns=columns)
	return (df)
