#!/usr/bin/env python
#python3 get_abYmod2.py /home/lilian/sync_project/Redundant_PDBs.txt /home/lilian/sync_project/input_Abymod_seq/ /home/lilian/sync_project/abYmod_structures/
#nohup python3 get_abYmod2.py /home/lilian/sync_project/Redundant_PDBs.txt /home/lilian/sync_project/input_Abymod_seq/ /home/lilian/sync_project/abYmod_structures/ &>output_abymod.log &
import sys
import os
#import shlex
#import subprocess
import shutil
#arg 0=filename; 1=Redundant_PDB file; 2=directory of .seq files; 3=directory to save model.pdb files;

def pass_commands(redundant_file, seq_directory, model_directory):
    file= open(redundant_file,"r")
    dic= {}
    list=[]
    for line in file:
        line=line.replace('\n','')
        list=line.split(',')
        for i in list:
            i.replace(" ","")
        list2=[]
        actual_seq_name=list[0]
        for i in list:
            if len(i)>4:
                new_code=i.split("_")
                nr=str(int(new_code[1])-1)
                letters=new_code[0].lower().replace(" ","")
                PDB_code=letters+nr
            else:
                PDB_code=i.lower().replace(" ","")
            list2.append(PDB_code)
        dic.update({actual_seq_name : list2})
    file.close()
    for key in dic:
        #command=shlex.split("abymod -v=3 -k=2 -exclude {} {} > {}".format(str(dic.get(key)).strip('[]').replace("'",""), (seq_directory+actual_seq_name+".seq"), (model_directory+actual_seq_name+".pdb.model")))
        #enter= subprocess.Popen(command, shell=True)
        to_exclude=str(dic.get(key)).strip('[]').replace("'","").replace(" ", "")
        key2=str(key).replace("'","")
        #print("abymod -v=3 -exclude={} -k=2 {} > {}".format(to_exclude, os.path.join(seq_directory+key2+".seq"), os.path.join(model_directory+key2+".pdb.model")))
        os.system("abymod -v=4 -excl100 -exclude={} {} > {} 2>{}".format(to_exclude, os.path.join(seq_directory,key2+".seq"), os.path.join(model_directory,key2+".pdb.model"), os.path.join(model_directory,key2+".log")))


def save_templates_seperately(seq_directory, model_directory):
    if os.path.exists(os.path.join(model_directory,"templates")):
        for filename in os.listdir(seq_directory):
            if ".tpl" in filename:
                source=os.path.join(seq_directory,filename)
                destination=os.path.join(model_directory,"templates",filename)
                shutil.move(source,destination)
    else:
        tpl_directory=os.mkdir(os.path.join(model_directory,"templates"))
        for filename in os.listdir(seq_directory):
            if ".tpl" in filename:
                source=os.path.join(seq_directory,filename)
                destination=os.path.join(model_directory,"templates", filename)
                shutil.move(source, destination)

pass_commands(sys.argv[1],sys.argv[2],sys.argv[3])
save_templates_seperately(sys.argv[2],sys.argv[3])
