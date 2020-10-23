#!/usr/bin/env python
file= open("/home/lilian/Desktop/structures/Redundant_PDBs.txt","r")
dic= {}
list=[]
for line in file:
    line=line.replace(',','')
    line=line.replace('\n','')
    list=line.split(' ')
    NR_PDB=list[0]
    del list[0]
    dic.update({NR_PDB : list})
file.close()
f = open("/home/lilian/sync_project/abYmod_commands.txt", "a")
for key in dic:
    string="abymod -v=3 -k=2 -exclude{} -autoexclude {} >{}".format(dic.get(key), (key+".pdb"), (key+"_model.pdb"))
    f.write(string + "\n")
f.close()

#then run in termianl using:
# bash home/lilian/sync_project/abYmod_commands.txt
#(maybe alter to save *_model in different directory)

#mv *.tpl /new_folder/*.tpl    =>save all template files in seperate folder
