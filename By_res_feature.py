#!/usr/bin/env python
import sys
import os
import numpy as np
import csv
#python3 By_res_feature.py ~/sync_project/abYmod_structures/ProFit_results/ ~/sync_project/Feature/ ~/sync_project/input_Abymod_seq/

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

def sort_by_res (profit_directory, feature_directory):
	fails=open(os.path.join(feature_directory,"fails"),'w+')
	for file in os.listdir(profit_directory):
		if os.stat(os.path.join(profit_directory,file)).st_size==0:
			filename=os.path.basename(file)
			fails.write(filename+"\n")
			continue
		write_line=[]
		with open(os.path.join(profit_directory,file),'r') as f:
			global local_AA
			global local_CA
			global global_AA
			global global_CA
			global name
			local_AA=[]
			local_CA=[]
			global_AA=[]
			global_CA=[]
			counter=0
			name=file.replace(".pdb.model.profit_by_res","")
			#print(name)
			for line in f:
				if ":" in line:
					if counter==0:
						local_AA.append(line)
						#print(line)
					if counter==1:
						local_CA.append(line)
					if counter==2:
						global_AA.append(line)
					if counter==3:
						global_CA.append(line)
				else:
					counter+=1
				worst_res(feature_directory)
    #return (local_AA, local_CA, global_AA, global_CA)
            #now you have the full lines for the file sorted-------------
def worst_res(feature_directory):
	separated=[]
	separated2=[]
	dics=[local_AA, local_CA, global_AA, global_CA]
	dics2=["local_AA", "local_CA", "global_AA", "global_CA"]
	for y in range (len(dics)):
		i=dics[y]
		for x in i:
			line=x.replace("\n","")
			line.strip("    ")
			separated.append(x)
		separated2 = [line.split() for line in separated]            #print(local_AA3)
		RMS={}
		for a in separated2:
			#print("separated",a)
			#print (name)
			rms=a[6]
			pos=a[3]
			pos=pos[1:]
			aa=a[1]
			RMS.update({rms:[pos,aa,dics2[y]]})
		if list(RMS.keys())==[]:
			print("empty:", name)
		worst_rms=max(list(RMS.keys()))
		write=[]
		'''binary classification of good and bad models
		if float(worst_rms) > 5:
			write.append(1)
		else:
			write.append(0)'''
		write.append(RMS[worst_rms][0])
		write.append(one_letter_code(RMS[worst_rms][1]))
		write.append(RMS[worst_rms][2])
		write.append(name)

		columns=["Pos","AA","RMS_type","ID"]
		if os.path.exists(os.path.join(feature_directory,"RMS_by_res_feature"+ ".csv")): #use feature_directory,"RMS_by_res_feature_csv_"+dics2[y] + ".csv") for separate files
			with open(os.path.join(feature_directory,"RMS_by_res_feature" + ".csv"), "a") as f:
				writer = csv.writer(f, delimiter=',')
				writer.writerow(write)
		else:
			with open(os.path.join(feature_directory,"RMS_by_res_feature"+ ".csv"),"w") as f:
				writer = csv.writer(f, delimiter=',')
				writer.writerow(columns)
				writer.writerow(write)


def Hydrophathy(feature_directory, input_seq_directory):
    Hydrophathy_index = {'A': 1.8, 'R': -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
                            "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2, "X": -0.5}#-0.5 is average
    columns=["Hydrophathy","Pos", "AA","ID"]
    with open(os.path.join(feature_directory,"hydropathy_index_feature_byres" + ".csv"),"w") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(columns)
    for filename in os.listdir(input_seq_directory):
        name= filename.replace(".seq","")
        with open(os.path.join(input_seq_directory,filename)) as infile:
            copy = False
            for line in infile:
                if "H95" in line:
                    copy = True
                if copy==True:
                    line2=line.split()
                    aminoacid=line2[1]
                    pos=line2[0][1:]
                    hydro_values=[Hydrophathy_index[aminoacid],pos,aminoacid,name]
                    with open(os.path.join(feature_directory,"hydropathy_index_feature_byres" + ".csv"),"a") as f:
                        writer = csv.writer(f, delimiter='\t')
                        writer.writerow(hydro_values)
                if "H102" in line:
                    copy = False


def surface(feature_directory, input_seq_directory):
    Accessible_surface_area= {'A': 44.1, 'R': 152.9, "N": 80.8, "D": 76.3, "C": 56.4, "Q": 100.6, "E": 99.2, "G": 0, "H": 98.2, "I": 90.9, "L": 92.8, "K": 139.1, "M": 95.3, "F": 107.4, "P": 79.5,
                            "S": 57.5,"T": 73.4, "W": 143.4, "Y": 119.1, "V": 73, "X":89}#89 is rounded average
    columns=["Surface","Pos","AA","ID"]
    with open(os.path.join(feature_directory,"accessible_surface_feature_byres" + ".csv"),"w") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(columns)
    for filename in os.listdir(input_seq_directory):
        name= filename.replace(".seq","")
        surface_values=[]
        with open(os.path.join(input_seq_directory,filename)) as infile:
            copy = False
            for line in infile:
                if "H95" in line:
                    copy = True
                if copy==True:
                    line2=line.split()
                    aminoacid=line2[1]
                    pos=line2[0][1:]
                    write=[Accessible_surface_area[aminoacid], pos, aminoacid, name]
                    with open(os.path.join(feature_directory,"accessible_surface_feature_byres" + ".csv"),"a") as f:
                        writer = csv.writer(f, delimiter='\t')
                        writer.writerow(write)
                if "H102" in line:
                    copy = False

sort_by_res(sys.argv[1], sys.argv[2])
Hydrophathy(sys.argv[2], sys.argv[3])
surface(sys.argv[2], sys.argv[3])
