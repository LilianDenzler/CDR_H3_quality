#!/usr/bin/env python
#python3 get_abYmod.py /home/lilian/sync_project/Redundant_PDBs.txt /home/lilian/sync_project/input_Abymod_seq/ /home/lilian/sync_project/abYmod_structures/
import sys
import os
import shlex
import subprocess
#arg 0=filename; 1=Redundant_PDB file; 2=directory of .seq files; 3=directory to save model.pdb files

def pass_commands(redundant_file, seq_directory, model_directory):
    file= open(redundant_file,"r")
    dic= {}
    list=[]
    for line in file:
        line=line[:-1]
        line=line.replace('\n','')
        list=line.split(', ')
        list2=[]
        actual_seq_name=list[0]
        print("actual name:{}".format(actual_seq_name))
        for i in list:
            PDB_code=(i[:i.find("_")])
            list2.append(PDB_code)
        dic.update({actual_seq_name : list2})
    file.close()
    bla=""
    for key in dic:
        #command=shlex.split("abymod -v=3 -k=2 -exclude {} {} > {}".format(str(dic.get(key)).strip('[]').replace("'",""), (seq_directory+actual_seq_name+".seq"), (model_directory+actual_seq_name+".pdb.model")))
        #enter= subprocess.Popen(command, shell=True)
        to_exclude=str(dic.get(key)).strip('[]').replace("'","").replace(" ", "")
        key2=str(key).replace("'","")
        print("abymod -v=3 -exclude={} -k=2 {} > {}".format(to_exclude, os.path.join(seq_directory+key2+".seq"), os.path.join(model_directory+key2+".pdb.model")))
        #os.system("abymod -v=3 -k=2 -exclude {} {} > {}".format(to_exclude, seq_directory+key2+".seq", model_directory+key2+".pdb.model"))

        '''def parse_abYmod_output():
            subprocess = subprocess.Popen("profit", "-f", "RMS_calc.txt", structure, shell=True, stdout=subprocess.PIPE)
            return_code = subprocess.poll()
            if return_code is not None:
              subprocess_return = subprocess.stdout.read()'''


if sys.argv[2].endswith("/")== False:
    sys.argv[2]=sys.argv[2]+"/"
if sys.argv[3].endswith("/")== False:
    sys.argv[3]=sys.argv[3]+"/"
pass_commands(sys.argv[1],sys.argv[2],sys.argv[3])
