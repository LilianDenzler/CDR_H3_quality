#!/usr/bin/env python
import sys
import os
#import numpy as np
import csv
#import pandas as pd
def accessibility(feature_directory, actual_directory):
    progress=0
    '''Accessible_surface_area= {'A': 44.1, 'R': 152.9, "N": 80.8, "D": 76.3, "C": 56.4, "Q": 100.6, "E": 99.2, "G": 0, "H": 98.2, "I": 90.9, "L": 92.8, "K": 139.1, "M": 95.3, "F": 107.4, "P": 79.5,
                            "S": 57.5,"T": 73.4, "W": 143.4, "Y": 119.1, "V": 73, "X":89}#89 is rounded average
    columns=["Surface","ID"]
    with open(os.path.join(feature_directory,"accessible_surface" + ".csv"),"w") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
    for filename in os.listdir(input_seq_directory):
        name= filename.replace(".seq","")
        surface_values=[]
        with open(os.path.join(input_seq_directory,filename)) as infile:
            copy = False
            surface=0
            for line in infile:
                if "H95" in line:
                    copy = True
                if copy==True:
                    line2=line.split()
                    aminoacid=line2[1]
                    pos=line2[0][1:]
                    surface+=Accessible_surface_area[aminoacid]
                if "H102" in line:
                    copy = False
            write=[surface,name]         
            with open(os.path.join(feature_directory,"accessible_surface" + ".csv"),"a") as f:
                writer = csv.writer(f)
                writer.writerow(write)'''
    columns=["Access","Relacc","Scacc","Screlacc","Access_avg","Relacc_avg","Scacc_avg","Screlacc_avg","ID"]
    with open(os.path.join(feature_directory,"accessibility" + ".csv"),"w") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
    for filename in os.listdir(actual_directory):
        nr_sad=0
        name= filename.replace(".pdb","")
        command=os.popen("pdbsolv -r stdout {}".format(os.path.join(actual_directory,filename))).readlines()# RESIDUE  AA   ACCESS  RELACC  SCACC   SCRELACC
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
            with open(os.path.join(feature_directory,"accessibility" + ".csv"),"a") as f:
                writer = csv.writer(f)
                writer.writerow(write)

        else:
            access_avg=access/counter
            relacc_avg=relacc/counter
            scacc_avg=scacc/counter
            screlacc_avg=screlacc/counter

            write=[access,relacc,scacc,screlacc, access_avg,relacc_avg, scacc_avg, screlacc_avg, name]
            with open(os.path.join(feature_directory,"accessibility" + ".csv"),"a") as f:
                writer = csv.writer(f)
                writer.writerow(write)
        progress+=1
        print(progress)
    f.close()

if __name__ == '__main__':
    accessibility(sys.argv[1], sys.argv[2])