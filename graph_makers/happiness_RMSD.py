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
#GRAPH: Sequence Similarity/Identity vs RMSD
#DO FITTING!!!!!!!!!!!!!!!!
#################################################################

data = pd.read_csv(sys.argv[1])
happy=data.Happiness_mean.tolist()
nr_sad=data.Nr_sad.tolist()
x=data.local_CA.tolist()
length=data.length.tolist()


del happy[0]
del nr_sad[0]
del x[0]
happy=[float(i) for i in happy]
nr_sad=[float(i) for i in nr_sad]
length=[int(i) for i in length]
nr_sad_norm=[]
for sad,l in zip(nr_sad, length):
	nr_sad_norm.append(sad*l)


x=[float(i) for i in x]
df1 = pd.DataFrame(list(zip(x, happy)), columns =['RMSD Local CA (Å)', 'mean happiness score']) 
df2 = pd.DataFrame(list(zip(x, nr_sad)), columns =['RMSD Local CA (Å)', "# residues with happiness score <0.5"]) 
df3 =pd.DataFrame(list(zip(x, nr_sad_norm)), columns =['RMSD Local CA (Å)', "# residues with happiness score <0.5*log(length)"]) 
ax=sns.jointplot(data=df1, x='mean happiness score', y='RMSD Local CA (Å)', kind='kde',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}})

plt.savefig("happy_RMSD.png")
plt.show()
sns.violinplot(data=df2, x="# residues with happiness score <0.5", y="RMSD Local CA (Å)",stat_func=stats.pearsonr, color="skyblue")
#plt.xlim([0,9])
plt.ylim([0,9])
plt.savefig("happy_nr_RMSD.png")
plt.show()

ax=sns.jointplot(data=df3, x='# residues with happiness score <0.5*log(length)', y='RMSD Local CA (Å)', kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}})

plt.savefig("happy_RMSD_reverseNorm.png")
plt.show()
 