#!/usr/bin/env python
import sys
import os
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
#GRAPH: Sequence Similarity/Identity vs RMSD
#DO FITTING!!!!!!!!!!!!!!!!
#################################################################
columns=["Happiness_mean", "Nr_sad","ID"]
data1 = pd.read_csv(sys.argv[1], names=columns)
columns2=["local_AA", "local_CA", "global_AA", "global_CA", "ID"]
data2 = pd.read_csv(sys.argv[2], names=columns2)
data=data1.merge(data2, on="ID")
happy=data.Happiness_mean.tolist()
nr_sad=data.Nr_sad.tolist()
x=data.local_CA.tolist()

del happy[0]
del nr_sad[0]
del x[0]
happy=[float(i) for i in happy]
nr_sad=[float(i) for i in nr_sad]
x=[float(i) for i in x]
df1 = pd.DataFrame(list(zip(x, happy)), columns =['RMSD Local CA (Å)', 'mean happiness score']) 
df2 = pd.DataFrame(list(zip(x, nr_sad)), columns =['RMSD Local CA (Å)', "# residues with mean happiness score <0.5"]) 
ax=sns.jointplot(data=df1, x='mean happiness score', y='RMSD Local CA (Å)', kind='kde',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}})

plt.savefig("happy_RMSD.png")
plt.show()
sns.violinplot(data=df2, x="# residues with mean happiness score <0.5", y="RMSD Local CA (Å)",stat_func=stats.pearsonr, color="skyblue")
#plt.xlim([0,9])
plt.ylim([0,9])
plt.savefig("happy_nr_RMSD.png")
plt.show()
 