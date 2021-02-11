#!/usr/bin/env python
#python3 get_abYmod2.py /home/lilian/sync_project/Redundant_PDBs.txt /home/lilian/sync_project/input_Abymod_seq/ /home/lilian/sync_project/abYmod_structures/
#nohup python3 get_abYmod2.py /home/lilian/sync_project/Redundant_PDBs.txt /home/lilian/sync_project/input_Abymod_seq/ /home/lilian/sync_project/abYmod_structures/ &>output_abymod.log &
import sys
import os
#import shlex
#import subprocess
import shutil


def pass_commands(seq_csv, upload_dir, ID):
	os.system("abymod -v=4 -excl100 {} > {} 2>{}".format(seq_csv, os.path.join(upload_dir,ID+".pdb.model"), os.path.join(upload_dir,ID+".log")))
	print("abymod -v=4 -excl100 {} > {} 2>{}".format(seq_csv, os.path.join(upload_dir,ID+".pdb.model"), os.path.join(upload_dir,ID+".log")))
	return (os.path.join(upload_dir,ID+".pdb.model"),os.path.join(upload_dir,ID+".log")) 
