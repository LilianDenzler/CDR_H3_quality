#!/usr/bin/env python
import sys
import os
import numpy as np
import csv
#sys.argv[1]=input_seq_directory,
def loop_redundancy(input_seq_directory):
    list_seq={}
    doubles={}
    for file in os.listdir(input_seq_directory:
        with open(os.path.join(input_seq_directory,file)) as f:
            name=file.replace(".seq","")
            copy = False
            sequence=""
            for line in f:
                if "H95" in line:
                    copy = True
                if copy==True:
                    line2=line.split()
                    aminoacid=line2[1]
                    sequence+=str(aminoacid)
                if "H102" in line:
                    copy = False
            if sequence in list_seq:
                doubles.update({sequence:name})
            else:
                list_seq.update({sequence:name})
        f.close()
    for keys,values in doubles.items():
        doubles[keys]=[doubles[keys],list_seq[keys]]




#just for checking:
'''for i in doubles.keys():
    for a in doubles[i]:
        with open(os.path.join(input_seq_directory,a+".seq")) as f:
            name=file.replace(".seq","")
            copy = False
            sequence=""
            for line in f:
                if "H95" in line:
                    copy = True
                if copy==True:
                    line2=line.split()
                    aminoacid=line2[1]
                    sequence+=str(aminoacid)
                if "H102" in line:
                    copy = False
            print(sequence)
        f.close()
'''
