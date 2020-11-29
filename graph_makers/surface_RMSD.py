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
columns=["Access", "Relacc", "Scacc", "Screlacc", "Access_avg", "Relacc_avg", "Scacc_avg", "Screlacc_avg","ID"]
data1 = pd.read_csv(sys.argv[1], names=columns)
columns2=["local_AA", "local_CA", "global_AA", "global_CA", "ID"]
data2 = pd.read_csv(sys.argv[2], names=columns2)
data=data1.merge(data2, on="ID")
access=data.Access.tolist()
relacc=data.Relacc.tolist()
scacc=data.Scacc.tolist()
screlacc=data.Screlacc.tolist()
access_avg=data.Access_avg.tolist()
relacc_avg=data.Relacc_avg.tolist()
scacc_avg=data.Scacc_avg.tolist()
screlacc_avg=data.Screlacc_avg.tolist()
x=data.local_CA.tolist()

del access[0]
del relacc[0]
del scacc[0]
del screlacc[0]
del access_avg[0]
del relacc_avg[0]
del scacc_avg[0]
del screlacc_avg[0]
del x[0]
access=[float(i) for i in access]
relacc=[float(i) for i in relacc]
scacc=[float(i) for i in scacc]
screlacc=[float(i) for i in screlacc]
access_avg=[float(i) for i in access_avg]
relacc_avg=[float(i) for i in relacc_avg]
scacc_avg=[float(i) for i in scacc_avg]
screlacc_avg=[float(i) for i in screlacc_avg]
x=[float(i) for i in x]


df1 = pd.DataFrame(list(zip(x, access)), columns =['RMSD Local CA (Å)', 'accessibility']) 
df2 = pd.DataFrame(list(zip(x, relacc)), columns =['RMSD Local CA (Å)', "relative accessibility"]) 
df3 = pd.DataFrame(list(zip(x, scacc)), columns =['RMSD Local CA (Å)', 'side-chain accessibility']) 
df4 = pd.DataFrame(list(zip(x, screlacc)), columns =['RMSD Local CA (Å)', "relative side-chain accessibility"]) 
df5 = pd.DataFrame(list(zip(x, access_avg)), columns =['RMSD Local CA (Å)', 'average accessibility']) 
df6 = pd.DataFrame(list(zip(x, relacc_avg)), columns =['RMSD Local CA (Å)', "average relative accessibility"]) 
df7 = pd.DataFrame(list(zip(x, scacc_avg)), columns =['RMSD Local CA (Å)', 'average side-chain accessibility']) 
df8 = pd.DataFrame(list(zip(x, screlacc_avg)), columns =['RMSD Local CA (Å)', "average relative side-chain accessibility"]) 
sns.jointplot(data=df1, x='accessibility', y='RMSD Local CA (Å)', kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,2200])
plt.savefig("surface_RMSD1.png")
sns.jointplot(data=df2, x="relative accessibility", y="RMSD Local CA (Å)", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,1200])
plt.savefig("surface_RMSD2.png")
sns.jointplot(data=df3, x="side-chain accessibility", y="RMSD Local CA (Å)", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,1800])
plt.savefig("surface_RMSD3.png")
sns.jointplot(data=df4, x="relative side-chain accessibility", y="RMSD Local CA (Å)", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,1500])
plt.savefig("surface_RMSD4.png")

sns.jointplot(data=df5, x='average accessibility', y='RMSD Local CA (Å)', kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,120])
plt.savefig("surface_RMSD1_avg.png")
sns.jointplot(data=df6, x="average relative accessibility", y="RMSD Local CA (Å)", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,80])
plt.savefig("surface_RMSD2_avg.png")
sns.jointplot(data=df7, x="average side-chain accessibility", y="RMSD Local CA (Å)", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,80])
plt.savefig("surface_RMSD3_avg.png")
sns.jointplot(data=df8, x="average relative side-chain accessibility", y="RMSD Local CA (Å)", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,65])
plt.savefig("surface_RMSD4_avg.png")
plt.show()