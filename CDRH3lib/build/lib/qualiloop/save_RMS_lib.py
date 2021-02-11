#!/usr/bin/env python
#python3 save_RMS ~/git/CDR_H3_quality/RMS_calc ~/sync_project/actual_PDBs_NR ~/sync_project/abYmod_structures ~/sync_project/Feature

import sys
import os
import csv

def pass_Profit_commands (script_name, actual_directory, model_directory, feature_directory):
    check_empty_model=os.popen("find {} -type f -empty".format(model_directory)).readlines()
    fails=open(os.path.join(actual_directory,"models"+".fails"),'w+')
    fails_list=[]
    data=[]
    for i in check_empty_model:
        i=os.path.basename(i)
        failed_model=i.replace(".pdb.model\n","")
        fails_list.append(failed_model)
        fails.write(i)

    for filename in os.listdir(model_directory):
        name=filename.replace(".pdb.model","")
        name=name.replace(" ","")
        print(name)
        if name in fails_list:
            print("yea")
        elif name =="ProFit_results":
            print("ok")
        else:
            command=os.popen("profit -f {} -h {} {}".format(script_name, os.path.join(actual_directory,name+".pdb"), os.path.join(model_directory,name+".pdb"+".model"))).readlines()
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
                    counter=1
                elif counter==1:
                    counter=0

            if len(list_write)>7:
                del list_write[2]
                del list_write[3]
                del list_write[3]
                del list_write[4]

            list_write.append(name)
            data.append(list_write)
    columns=["local_AA","local_CA","global_AA","global_CA","ID"]
    df = pd.DataFrame(data, columns=columns)
    print("get_loop_seq:",df.shape)
    return (df)

def run(actual_directory, model_directory, feature_directory):
    pass_Profit_commands("./RMS_calc", actual_directory, model_directory, feature_directory)
