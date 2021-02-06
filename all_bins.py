#!/usr/bin/env python
import sys
import os
import pandas as pd
import math
#import save_RMS_lib
#import get_abYmod2_lib
#import seq_parser_lib
from pathlib import Path
import numpy as np
import epitopepredict as ep
from skimage.feature import hog
from skimage import data, exposure
from sklearn.decomposition import TruncatedSVD
import collections
from sklearn import manifold, datasets
from functools import partial

def RMSD_binary(actual_directory, model_directory, feature_directory, threshold,feature_csv):
    columns=["local_AA_bin"+str(threshold),"local_CA_bin"+str(threshold),"global_AA_bin"+str(threshold), "global_CA_bin"+str(threshold), "ID"]
    data=[]
    RMSD_path=Path(feature_csv)
    if os.path.exists(RMSD_path) ==True:
        df=pd.read_csv(RMSD_path, header=0)
    local_AA= df.local_AA.tolist()
    local_CA= df.local_CA.tolist()
    global_AA= df.global_AA.tolist()
    global_CA=df.global_CA.tolist()
    name=df.ID.tolist()
    for (a,b,c,d,n) in zip(local_AA, local_CA, global_AA, global_CA, name):
        a=float(a)
        b=float(b)
        c=float(c)
        d=float(d)
        if a>=threshold:
            a=1
        else:
            a=0
        if b>=threshold:
            b=1
        else:
            b=0
        if c>=threshold:
            c=1
        else:
            c=0
        if d>=threshold:
            d=1
        else:
            d=0
        write=[a,b,c,d,n]
        data.append(write)
    dfnew = pd.DataFrame(data, columns=columns)
    print(dfnew)
    return (dfnew)

def RMSD_nom(actual_directory, model_directory, feature_directory, RMSD_file,mode,stop):
    columns=["local_AA_nom","local_CA_nom","global_AA_nom", "global_CA_nom", "ID"]
    data=[]
    RMSD_path=Path(RMSD_file)
    if os.path.exists(RMSD_path) ==True:
        df=pd.read_csv(RMSD_path, header=0)
    else:
        input("AAAAH")
    local_AA= df.local_AA.tolist()
    local_CA= df.local_CA.tolist()
    global_AA= df.global_AA.tolist()
    global_CA=df.global_CA.tolist()
    name=df.ID.tolist()
    if mode =="half":
    	threshlow=list(np.linspace(0,4,9))
    	threshup=list(np.linspace(0.5,4,8))
    	print("threshlow",threshlow)
    	print("thresholdup",threshup)
    	values=[">"+str(x-0.5)+"<"+str(x) for x in threshup]
    	values.append(">"+str(4.5))
    if mode=="full":
    	threshlow=list(range(0,int(stop),2))
    	threshup=list(range(2,int(stop),2))
    	print("threshlow",threshlow)
    	print("thresholdup",threshup)
    	values=[">"+str(x-2)+"<"+str(x) for x in threshup]
    	values.append(">"+str(int(stop)-2))
    threshup.append(10000)
    print(values)
    for (a,b,c,d,n) in zip(local_AA, local_CA, global_AA, global_CA, name):
        for (threshold_low,threshold_up,val) in list(zip(threshlow,threshup,values)):
            if type(a)==str:
                pass
            elif a>=threshold_low and a<threshold_up:
                a=val
            if type(b)==str:
                pass
            elif b>=threshold_low and b<threshold_up:
                b=val
            if type(c)==str:
                pass
            elif c>=threshold_low and c<threshold_up:
                c=val
            if type(d)==str:
                pass
            elif d>=threshold_low and d<threshold_up:
                d=val

        write=[a,b,c,d,n]
        data.append(write)
    dfnew = pd.DataFrame(data, columns=columns)
    print("RMSD_nom:",dfnew.shape)
    dfnew.to_csv(os.path.join(feature_directory,"bins"+mode+".csv"), index=False)
    print(dfnew)
    return (dfnew)

def merge(actual_directory, model_directory, feature_directory, feature_csv):
	mode=input("for half-bins (i.e. 0-0.5, 0.5-1,..>4.5 type half, if you want full bins (i.e 0-1, 1-2,...>5) type full")
	if mode=="half":
		threshlow=list(np.linspace(0,4,9))
	if mode=="full":
		stop=input("number of bins")
		threshlow=list(range(2,int(stop)+2,2))
	bin_dict={}
	a=pd.read_csv(feature_csv)
	#a=a.drop(['local_AA_bin',"global_AA_bin", "local_CA_bin","global_CA_bin", 'local_AA_nom',"global_AA_nom", "local_CA_nom","global_CA_nom"],axis=1)
	print(a.shape)
	for i in threshlow:
		df=RMSD_binary(actual_directory, model_directory, feature_directory, i,feature_csv)
		bin_dict.update({i:df})
	df_nom=RMSD_nom(actual_directory, model_directory, feature_directory, feature_csv,mode,stop)
	bin_dict.update({"nom":df_nom})
	for i in bin_dict.keys():
		b=bin_dict[i]
		merged = a.concat(b, on='ID')
		print("merged:", merged.shape)
		a=merged
	a.to_csv(os.path.join(feature_directory,"basics"+".csv"), index=False)
	print(a)

merge(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
#python3 all_bins.py ~/sync_project/actual_PDBs_NR/ ~/sync_project/abymod_structures2/ ~/sync_project/Feature/ ~/sync_project/Feature/new_merged.csv