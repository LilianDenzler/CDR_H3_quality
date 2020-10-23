#!/usr/bin/env python
import sys
import os
#arg 0=filename; 1=Redundant_PDB file; 2=directory of .seq files; 3=directory to save model.pdb files
def pass_commands(redundant_file, seq_directory, model_directory):
    file= open(redundant_file,"r")
    dic= {}
    list=[]
    for line in file:
        line=line.replace(',','')
        line=line.replace('\n','')
        list=line.split(' ')
        NR_PDB=list[0]
        dic.update({NR_PDB : list})
    file.close()
    bla=""
    for key in dic:
        os.system("abymod -v=3 -k=2 -exclude{} {} > {}".format(str(dic.get(key)).strip('[]'), (seq_directory +"/"+key+".pdb"), (model_directory+"/"+key+".model.pdb")))

pass_commands(sys.argv[1],sys.argv[2],sys.argv[3])
