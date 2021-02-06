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
data = pd.read_csv(sys.argv[1], header=0)
access=data.Access.tolist()
relacc=data.Relacc.tolist()
scacc=data.Scacc.tolist()
screlacc=data.Screlacc.tolist()
access_avg=data.Access_avg.tolist()
relacc_avg=data.Relacc_avg.tolist()
scacc_avg=data.Scacc_avg.tolist()
screlacc_avg=data.Screlacc_avg.tolist()
x=data.local_CA.tolist()
length=data.length.tolist()
length=[int(x) for x in length]

#sns.set()
sns.jointplot(data=data, x='Access_avg', y='local_CA', color=None,kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}},scatter_kws={'alpha':0.0})
ax=sns.scatterplot(data=data, x='Access_avg', y='local_CA',hue="length", palette='viridis',alpha=0.5, s=50)
sm = plt.cm.ScalarMappable(cmap="viridis")
sm.set_array([])
plt.ylim([0,9])
plt.xlim([0,100])
ax.get_legend().remove()
#ax.figure.colorbar(sm)
plt.savefig("surface_Acc_avg.png")

sns.jointplot(data=data, x='Relacc_avg', y='local_CA', color=None,kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}},scatter_kws={'alpha':0.0})
ax=sns.scatterplot(data=data, x='Relacc_avg', y='local_CA',hue="length", palette='viridis',alpha=0.5, s=50)
sm = plt.cm.ScalarMappable(cmap="viridis")
sm.set_array([])
plt.ylim([0,9])
plt.xlim([0,80])
ax.get_legend().remove()
ticks=[0,2,4,6,8,10, 12,14,16,18,20,22,24,26,28]
nticks=[((x - min(length)) / (max(length) - min(length))) for x in ticks]
ticks=[str(x) for x in ticks]
#cbar=ax.figure.colorbar(sm, ticks=nticks, orientation='horizontal')
#cbar.ax.set_xticklabels(ticks) 
plt.savefig("surface_Relacc_avg.png")

sns.jointplot(data=data, x='Scacc_avg', y='local_CA', color=None,kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}},scatter_kws={'alpha':0.0})
ax=sns.scatterplot(data=data, x='Scacc_avg', y='local_CA',hue="length", palette='viridis',alpha=0.5, s=50)
sm = plt.cm.ScalarMappable(cmap="viridis")
sm.set_array([])
plt.ylim([0,9])
plt.xlim([0,80])
ax.get_legend().remove()
ticks=[0,2,4,6,8,10, 12,14,16,18,20,22,24,26,28]
nticks=[((x - min(length)) / (max(length) - min(length))) for x in ticks]
ticks=[str(x) for x in ticks]
#cbar=ax.figure.colorbar(sm, ticks=nticks, orientation='horizontal')
#cbar.ax.set_xticklabels(ticks) 
plt.savefig("surface_Scacc_avg.png")

sns.jointplot(data=data, x='Screlacc_avg', y='local_CA', color=None,kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}},scatter_kws={'alpha':0.0})
ax=sns.scatterplot(data=data, x='Screlacc_avg', y='local_CA',hue="length", palette='viridis',alpha=0.5, s=50)
sm = plt.cm.ScalarMappable(cmap="viridis")
sm.set_array([])
plt.ylim([0,9])
plt.xlim([0,65])
ax.get_legend().remove()
#ax.figure.colorbar(sm)
plt.savefig("surface_Screlacc_avg.png")

sns.jointplot(data=data, x='Access', y='local_CA', color=None,kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}},scatter_kws={'alpha':0.0})
ax=sns.scatterplot(data=data, x='Access', y='local_CA',hue="length", palette='viridis',alpha=0.5, s=50)
sm = plt.cm.ScalarMappable(cmap="viridis")
sm.set_array([])
plt.ylim([0,9])
plt.xlim([0,1200])
ax.get_legend().remove()
#ax.figure.colorbar(sm)
plt.savefig("surface_Acc.png")

sns.jointplot(data=data, x='Relacc', y='local_CA', color=None,kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}},scatter_kws={'alpha':0.0})
ax=sns.scatterplot(data=data, x='Relacc', y='local_CA',hue="length", palette='viridis',alpha=0.5, s=50)
sm = plt.cm.ScalarMappable(cmap="viridis")
sm.set_array([])
plt.ylim([0,9])
plt.xlim([0,1200])
ax.get_legend().remove()
ticks=[0,2,4,6,8,10, 12,14,16,18,20,22,24,26,28]
nticks=[((x - min(length)) / (max(length) - min(length))) for x in ticks]
ticks=[str(x) for x in ticks]
#cbar=ax.figure.colorbar(sm, ticks=nticks, orientation='horizontal')
#cbar.ax.set_xticklabels(ticks) 
plt.savefig("surface_Relacc.png")

sns.jointplot(data=data, x='Scacc', y='local_CA', color=None,kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}},scatter_kws={'alpha':0.0})
ax=sns.scatterplot(data=data, x='Scacc', y='local_CA',hue="length", palette='viridis',alpha=0.5, s=50)
sm = plt.cm.ScalarMappable(cmap="viridis")
sm.set_array([])
plt.ylim([0,9])
plt.xlim([0,1800])
ax.get_legend().remove()
ticks=[0,2,4,6,8,10, 12,14,16,18,20,22,24,26,28]
nticks=[((x - min(length)) / (max(length) - min(length))) for x in ticks]
ticks=[str(x) for x in ticks]
#cbar=ax.figure.colorbar(sm, ticks=nticks, orientation='horizontal')
#cbar.ax.set_xticklabels(ticks) 
plt.savefig("surface_Scacc.png")

sns.jointplot(data=data, x='Screlacc', y='local_CA', color=None,kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}},scatter_kws={'alpha':0.0})
ax=sns.scatterplot(data=data, x='Screlacc', y='local_CA',hue="length", palette='viridis',alpha=0.5, s=50)
sm = plt.cm.ScalarMappable(cmap="viridis")
sm.set_array([])
plt.ylim([0,9])
plt.xlim([0,1500])
ax.get_legend().remove()
#ax.figure.colorbar(sm)
plt.savefig("surface_Screlacc.png")


'''



sns.jointplot(data=data, x="Access", y="local_CA", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,1200])
plt.savefig("surface_RMSD1.png")

sns.jointplot(data=data, x="Relacc", y="local_CA", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,1200])
plt.savefig("surface_RMSD2.png")
sns.jointplot(data=data, x="Scacc", y="local_CA", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,1800])
plt.savefig("surface_RMSD3.png")
sns.jointplot(data=data, x="Screlacc", y="local_CA", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,1500])
plt.savefig("surface_RMSD4.png")

sns.jointplot(data=data, x='Access_avg', y='local_CA', kind='reg',stat_func=stats.pearsonr, joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,120])
plt.savefig("surface_RMSD1_avg.png")
sns.jointplot(data=data, x="Relacc_avg", y="local_CA", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,80])
plt.savefig("surface_RMSD2_avg.png")
sns.jointplot(data=data, x="Scacc_avg", y="local_CA", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,80])
plt.savefig("surface_RMSD3_avg.png")
sns.jointplot(data=data, x="Screlacc_avg", y="local_CA", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.ylim([0,9])
plt.xlim([0,65])
plt.savefig("surface_RMSD4_avg.png")
plt.show()'''