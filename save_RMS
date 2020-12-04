#!/usr/bin/env python
#python3 save_RMS ~/git/CDR_H3_quality/RMS_calc ~/sync_project/actual_PDBs_NR ~/sync_project/abYmod_structures ~/sync_project/Feature

import sys
import os
import csv

def pass_Profit_commands (script_name, actual_directory, model_directory, feature_directory):
    check_empty_model=os.popen("find {} -type f -empty".format(model_directory)).readlines()
    fails=open(os.path.join(actual_directory,"fails"),'w+')
    fails_list=[]
    for i in check_empty_model:
        i=os.path.basename(i)
        fails.write(i)
        failed_model=i.replace(".pdb.model\n","")
        fails_list.append(failed_model)

    for filename in os.listdir(model_directory):
        name=filename.replace(".pdb.model","")
        name=name.replace(" ","")
        print(name)
        if name in fails_list:
            print("yea")
        elif name =="ProFit_results":
            print("ok")
        else:
            command=os.popen("./profit -f {} -h {} {}".format(script_name, os.path.join(actual_directory,name+".pdb"), os.path.join(model_directory,name+".pdb"+".model"))).readlines()
            #RMS=open(os.path.join(model_directory,"ProFit_results",filename + ".profit"),'w+')
            if not os.path.exists(os.path.join(feature_directory,"ProFit_results",filename+ ".profit_by_res")):
                os.makedirs(os.path.join(feature_directory,"ProFit_results",filename+ ".profit_by_res"))
            By_res=open(os.path.join(feature_directory,"ProFit_results",filename+ ".profit_by_res"),'w+')

            counter=0
            list_write=[]
            for line in command:
                print(line)
                columns=line.split("    ")
                if len(columns)< 2 and 'RMS:' in line:
                    line=line.replace(" ","")
                    line=line.replace("\n","")
                    list_write.append(line[4:])
                if len(columns)>=3 and 'RMS:' in line:
                    By_res.write(line)
                    counter=1
                elif counter==1:
                    By_res.write(" \n")
                    counter=0

            if len(list_write)>7:
                del list_write[2]
                del list_write[3]
                del list_write[3]
                del list_write[4]

            list_write.append(name)
            print(list_write)
            #print(list_write)
            columns=["local_AA","local_CA","global_AA","global_CA","ID"]
            if os.path.exists(os.path.join(feature_directory,"RMS_feature_csv" + ".csv")):
                with open(os.path.join(feature_directory,"RMS_feature_csv" + ".csv"), "a") as f:
                    writer = csv.writer(f)
                    writer.writerow(list_write)
            else:
                with open(os.path.join(feature_directory,"RMS_feature_csv" + ".csv"),"w") as f:
                    writer = csv.writer(f)
                    writer.writerow(columns)
                    writer.writerow(list_write)


    By_res.close()
    f.close()


'''
The .profit output file will have the following lines:
RMS local atom-to-ATOM
RMS local Ca-to-Ca
RMS global atom-to-atom just loop
RMS global Ca-to-Ca just loop
PDB ID

The .profit_by_res file will have the following lines:
RMS local atom-to-ATOM by residue
RMS local Ca-to-Ca by residue
RMS global atom-to-atom by residue
RMS global Ca-to-Ca by residue
'''

pass_Profit_commands(sys.argv[1],sys.argv[2],sys.argv[3], sys.argv[4])
