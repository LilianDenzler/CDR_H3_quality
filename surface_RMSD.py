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

#################################################################
#GRAPH: Sequence Similarity/Identity vs RMSD
#DO FITTING!!!!!!!!!!!!!!!!
#################################################################
columns=["Access", "Relacc", "Scacc", "Screlacc", "ID"]
data1 = pd.read_csv(sys.argv[1], names=columns)
columns2=["local_AA", "local_CA", "global_AA", "global_CA", "ID"]
data2 = pd.read_csv(sys.argv[2], names=columns2)
data=data1.merge(data2, on="ID")
access=data.Access.tolist()
relacc=data.Relacc.tolist()
scacc=data.Scacc.tolist()
screlacc=data.Screlacc.tolist()
x=data.local_CA.tolist()

del access[0]
del relacc[0]
del scacc[0]
del screlacc[0]
del x[0]
access=[float(i) for i in access]
relacc=[float(i) for i in relacc]
scacc=[float(i) for i in scacc]
screlacc=[float(i) for i in screlacc]
x=[float(i) for i in x]

print(type(access[2]))

df1 = pd.DataFrame(list(zip(x, access)), columns =['RMSD Local CA (Å)', 'accessibility']) 
df2 = pd.DataFrame(list(zip(x, relacc)), columns =['RMSD Local CA (Å)', "relative accessibility"]) 
df3 = pd.DataFrame(list(zip(x, scacc)), columns =['RMSD Local CA (Å)', 'side-chain accessibility']) 
df4 = pd.DataFrame(list(zip(x, screlacc)), columns =['RMSD Local CA (Å)', "relative side-chain accessibility"]) 
sns.jointplot(data=df1, x='accessibility', y='RMSD Local CA (Å)', kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,260])
plt.savefig("surface_RMSD1.png")
sns.jointplot(data=df2, x="relative accessibility", y="RMSD Local CA (Å)", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,120])
plt.savefig("surface_RMSD2.png")
sns.jointplot(data=df3, x="side-chain accessibility", y="RMSD Local CA (Å)", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,200])
plt.savefig("surface_RMSD3.png")
sns.jointplot(data=df4, x="relative side-chain accessibility", y="RMSD Local CA (Å)", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,120])
plt.savefig("surface_RMSD4.png")
plt.show()