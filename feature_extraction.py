#!/usr/bin/env python
import sys
import os
import numpy as np
import csv

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

def get_loop_charge:
    dic = {'C': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','XAA': 'X', 'UNK':'X'}

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

def get_

#python3 feature_extraction.py ~/sync_project/Feature/ ~/sync_project/input_Abymod_seq/

if __name__=="__main__":
  get_loop_length(sys.argv[1],sys.argv[2])
