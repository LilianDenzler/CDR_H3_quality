#!/usr/bin/env python
import sys
import os
import math
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import optimize
import seaborn as sns
from sklearn import preprocessing

#################################################################
#GRAPH: Different encodings vs RMSD
#DO FITTING!!!!!!!!!!!!!!!!
#################################################################

data = pd.read_csv(sys.argv[1], header=0)
Mod_LLE1=data.Mod_LLE1.tolist()
Mod_LLE2=data.Mod_LLE2.tolist()
Mod_LLE3=data.Mod_LLE3.tolist()
Mod_LLE4=data.Mod_LLE4.tolist()
Mod_LLE5=data.Mod_LLE5.tolist()
Mod_LLE6=data.Mod_LLE6.tolist()
Isomap1=data.Isomap1.tolist()
Isomap2=data.Isomap2.tolist()
Isomap3=data.Isomap3.tolist()
Isomap4=data.Isomap4.tolist()
Isomap5=data.Isomap5.tolist()
Isomap6=data.Isomap6.tolist()
blosum62_val1=data.blosum62_val1.tolist()
blosum62_val2=data.blosum62_val2.tolist()
blosum62_val3=data.blosum62_val3.tolist()
blosum62_val4=data.blosum62_val4.tolist()
blosum62_val5=data.blosum62_val5.tolist()
blosum62_val6=data.blosum62_val6.tolist()
protrusion=data.protrusion.tolist()
y=data.local_CA.tolist()

#happy=[float(i) for i in happy]

#x=[float(i) for i in x]
fig = plt.figure()
plt.xlabel('Mod_LLE')
plt.ylabel('local_CA')
plt.plot(Mod_LLE1, y,'o',markersize=5)
plt.plot(Mod_LLE2, y,'ro',markersize=5)
plt.plot(Mod_LLE3, y,'bo',markersize=5)
plt.plot(Mod_LLE4, y,'yo',markersize=5)
plt.plot(Mod_LLE5, y,'go',markersize=5)
plt.plot(Mod_LLE6, y,'mo',markersize=5)
plt.show()

fig = plt.figure()
plt.xlabel('Isomap')
plt.ylabel('local_CA')
plt.plot(Isomap1, y,'o',markersize=5)
plt.plot(Isomap2, y,'ro',markersize=5)
plt.plot(Isomap3, y,'bo',markersize=5)
plt.plot(Isomap4, y,'yo',markersize=5)
plt.plot(Isomap5, y,'go',markersize=5)
plt.plot(Isomap6, y,'mo',markersize=5)
plt.show()

fig = plt.figure()
sns.jointplot(data=data, x='local_CA', y='protrusion', kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}})
plt.show()

fig = plt.figure()
plt.xlabel('blosum62')
plt.ylabel('local_CA')
plt.plot(blosum62_val1, y,'o',markersize=5)
plt.plot(blosum62_val2, y,'ro',markersize=5)
plt.plot(blosum62_val3, y,'bo',markersize=5)
plt.plot(blosum62_val4, y,'yo',markersize=5)
plt.plot(blosum62_val5, y,'go',markersize=5)
plt.plot(blosum62_val6, y,'mo',markersize=5)
plt.show()
################################################



fig = plt.figure()
sns.jointplot(data=data, x='local_CA', y='Mod_LLE1', kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}})
sns.jointplot(data=data, x='local_CA', y='Mod_LLE2', kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}})
sns.jointplot(data=data, x='local_CA', y='Mod_LLE3', kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}})
sns.jointplot(data=data, x='local_CA', y='Mod_LLE4', kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}})
sns.jointplot(data=data, x='local_CA', y='Mod_LLE5', kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}})
sns.jointplot(data=data, x='local_CA', y='Mod_LLE6', kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}})
#plt.savefig("Mod_LLE_RMSD.png")
plt.show()

#sns.violinplot(data=df2, x="# residues with happiness score <0.5", y="RMSD Local CA (Å)",stat_func=stats.pearsonr, color="skyblue")
#plt.xlim([0,9])
#plt.ylim([0,9])
#plt.savefig("happy_nr_RMSD.png")
#plt.show()

#ax=sns.jointplot(data=df3, x='# residues with happiness score <0.5*log(length)', y='RMSD Local CA (Å)', kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}})

#plt.savefig("happy_RMSD_reverseNorm.png")
#plt.show()
 