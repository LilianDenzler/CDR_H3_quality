#!/usr/bin/env python
import sys
import os
import numpy as np
import csv

def extract_seq(feature_directory,input_seq_directory):
    pass
def write_csv(feature_directory):
    pass


def get_loop_length (feature_directory,input_seq_directory):
    columns=["length","ID"]
    with open(os.path.join(feature_directory,"length_feature_csv" + ".csv"),"w") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(columns)
    for filename in os.listdir(input_seq_directory):
        name= filename.replace(".seq","")
        with open(os.path.join(input_seq_directory,filename)) as infile:
            copy = False
            count= 0
            for line in infile:
                if "H95" in line:
                    count+=1
                    copy = True
                    continue
                elif "H102" in line:
                    count+=1
                    copy = False
                    continue
                elif copy:
                    count+=1
        write=[str(count)]
        write.append(str(name))
        print(write)
        with open(os.path.join(feature_directory,"length_feature_csv" + ".csv"),"a") as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(write)

def get_loop_charge(feature_directory,input_seq_directory):
    dic = {"D":-1, "K": -1,"R": 1,'E': 1, 'H':1}
    columns=["charge","ID"]
    with open(os.path.join(feature_directory,"charge_feature_csv" + ".csv"),"w") as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(columns)
    for filename in os.listdir(input_seq_directory):
        name= filename.replace(".seq","")
        with open(os.path.join(input_seq_directory,filename)) as infile:
            copy = False
            aminoacid=[]
            charge=0
            for line in infile:
                if "H95" in line:
                    copy = True
                if copy==True:
                    line2=line.split()
                    #print(line2)
                    aminoacid=line2[1]
                    #print ("AA:",aminoacid)
                    if aminoacid in dic:
                        charge+=dic[aminoacid]
                        print(name, charge, aminoacid, line2[0])
                if "H102" in line:
                    copy = False
        write=[str(charge)]
        write.append(str(name))
        #print(write)
        with open(os.path.join(feature_directory,"charge_feature_csv" + ".csv"),"a") as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(write)

def hydrophobicity(feature_directory,input_seq_directory):
    Hydrophathy_index = {'A': 1.8, 'R': -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
                            "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2, "X": -0.5}#-0.5 is average
    columns=["Hydrophathy","Pos", "AA","ID"]
    with open(os.path.join(feature_directory,"hydropathy_index_feature" + ".csv"),"w") as f:
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
                    with open(os.path.join(feature_directory,"hydropathy_index_feature" + ".csv"),"a") as f:
                        writer = csv.writer(f, delimiter='\t')
                        writer.writerow(hydro_values)
                if "H102" in line:
                    copy = False


def surface(feature_directory, input_seq_directory):
    Accessible_surface_area= {'A': 44.1, 'R': 152.9, "N": 80.8, "D": 76.3, "C": 56.4, "Q": 100.6, "E": 99.2, "G": 0, "H": 98.2, "I": 90.9, "L": 92.8, "K": 139.1, "M": 95.3, "F": 107.4, "P": 79.5,
                            "S": 57.5,"T": 73.4, "W": 143.4, "Y": 119.1, "V": 73, "X":89}#89 is rounded average
    columns=["Surface","Pos","AA","ID"]
    with open(os.path.join(feature_directory,"accessible_surface_feature" + ".csv"),"w") as f:
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
                    with open(os.path.join(feature_directory,"accessible_surface_feature" + ".csv"),"a") as f:
                        writer = csv.writer(f, delimiter='\t')
                        writer.writerow(write)
                if "H102" in line:
                    copy = False


def Z-scale(feature_directory,input_seq_directory):
    pass
def VHSE_scale(feature_directory,input_seq_directory):
    pass
def template_accuracy(feature_directory,input_seq_directory):
    pass
def bulged_non-bulged(feature_directory,input_seq_directory):
    pass
def nr_contacts(feature_directory,input_seq_directory):
    pass
def seq_ransomisation(feature_directory,input_seq_directory):
    pass
    #also information entropy
def BLOSUM_62():
    pass
def BLOSUM_50():
    pass

#also, how many AA of same category after one another
'''
Charged (side chains often form salt bridges):
• Arginine - Arg - R
• Lysine - Lys - K
• Aspartic acid - Asp - D
• Glutamic acid - Glu - E

Polar (form hydrogen bonds as proton donors or acceptors):
• Glutamine - Gln - Q
• Asparagine - Asn - N
• Histidine - His - H
• Serine - Ser - S
• Threonine - Thr - T
• Tyrosine - Tyr - Y
• Cysteine - Cys - C

Amphipathic (often found at the surface of proteins or lipid membranes, sometimes also classified as polar):
• Tryptophan - Trp - W
• Tyrosine - Tyr - Y
• Methionine - Met - M (may function as a ligand to metal ions)

Hydrophobic (normally buried inside the protein core):
• Alanine - Ala - A
• Isoleucine - Ile - I
• Leucine - Leu - L
• Methionine - Met - M
• Phenylalanine - Phe - F
• Valine - Val - V
• Proline - Pro - P
• Glycine - Gly - G
'''

def merge_csv():
    import pandas as pd
    a = pd.read_csv("filea.csv")
    b = pd.read_csv("fileb.csv")
    b = b.dropna(axis=1)
    merged = a.merge(b, on='title')
    merged.to_csv("output.csv", index=False)
#python3 feature_extraction.py ~/sync_project/Feature/ ~/sync_project/input_Abymod_seq/

if __name__=="__main__":
  #get_loop_length(sys.argv[1],sys.argv[2])
  get_loop_charge(sys.argv[1],sys.argv[2])
