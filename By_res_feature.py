#!/usr/bin/env python
import sys
import os
import numpy as np
import csv

def pass_Profit_commands (profit_directory, feature_directory):
    for file in os.listdir(profit_directory):
        if os.stat(os.path.join(profit_directory,file)).st_size==0:
            continue
        write_line=[]
        with open(os.path.join(profit_directory,file),'r') as f:
            local_AA=[]
            local_CA=[]
            global_AA=[]
            global_CA=[]
            local_AA2=[]
            local_CA2=[]
            global_AA2=[]
            global_CA2=[]
            counter=0
            name=file.replace(".pdb.model.profit_by_res","")
            #print(name)
            for line in f:
                if ":" in line:
                    if counter==0:
                        local_AA.append(line)
                        print(line)
                    if counter==1:
                        local_CA.append(line)
                    if counter==2:
                        global_AA.append(line)
                    if counter==3:
                        global_CA.append(line)
                else:
                    counter+=1
            #now you have the full lines for the file sorted

            #print(local_AA)
            for i in local_AA:
                i=i.replace("\n","")
                i.strip("    ")
                #i.split('    ')
                local_AA2.append(i)
            local_AA3 = [i.split() for i in local_AA2]
            #print(local_AA3)
            local_AA_RMS={}
            for i in local_AA3:
                rms=i[6]
                pos=i[3]
                pos=pos[1:]
                aa=i[1]
                local_AA_RMS.update({rms:[pos,aa]})
            #print(local_AA_RMS)
            worst_rms=max(list(local_AA_RMS.keys()))
            #print("worst",worst_rms)
            write=[]
            if float(worst_rms) > 5:
                write.append(1)
            else:
                write.append(0)
            write.append(local_AA_RMS[worst_rms][0])
            write.append(local_AA_RMS[worst_rms][1])
            #print(write)

        columns=["RMSD","Pos","AA"]
        if os.path.exists(os.path.join(feature_directory,"RMS_by_res_feature_csv_localAA" + ".csv")):
            with open(os.path.join(feature_directory,"RMS_by_res_feature_csv_localAA" + ".csv"), "a") as f:
                writer = csv.writer(f, delimiter=',')
                writer.writerow(write)
        else:
            with open(os.path.join(feature_directory,"RMS_by_res_feature_csv_localAA" + ".csv"),"w") as f:
                writer = csv.writer(f, delimiter=',')
                writer.writerow(columns)
                writer.writerow(write)


pass_Profit_commands(sys.argv[1], sys.argv[2])
