#!/usr/bin/env python
import sys
import os
import pandas
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import optimize
import seaborn as sns
#################################################################
#GRAPH: CDRH3 Loop Ca-RMSD vs Loop Length
#################################################################
columns=["local_AA","local_CA","global_AA","global_CA","ID", "length", "identity", "similarity", "template", "target"]
data = pandas.read_csv(sys.argv[1], names=columns)
x=data.length.tolist()
a=data.local_CA.tolist()
b=data.global_CA.tolist()
del x[0]
del a[0]
del b[0]
x = [g for g in x if str(g) != 'nan']
x = [float(g) for g in x]
a = [g for g in a if str(g) != 'nan']
a = [float(g) for g in a]
b = [g for g in b if str(g) != 'nan']
b = [float(g) for g in b]

r1=stats.pearsonr(x, a) #	(Pearson’s correlation coefficient, 2-tailed p-value)
r2=stats.pearsonr(x, b)

m1, b1 = np.polyfit(x, a, 1)
m2, b2 = np.polyfit(x, b, 1)


fig, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(50, 50))# sharex=True, sharey=True
plt.subplots_adjust(wspace=0.5, hspace=0.5)
fig.suptitle('CDRH3 Loop RMSD vs Loop Length', size="20", fontweight='bold')
ax1.set_ylabel('Local C-alpha RMSD (Å)')
ax1.set_xlabel('Loop Length (#Residues)')
ax1.plot(x, a,'o',markersize=1)
ax1.text(0,6.5, "Pearson correlation coefficient:{}".format(round(r1[0],2)), fontsize=10, fontweight='bold')
ax1.plot(np.unique(x), m1*np.unique(x) + b1, label="y={}x+{}".format(round(m1, 2), round(b1, 2)))
ax1.legend(loc="center left", bbox_to_anchor=(1, 0.5))

ax2.set_ylabel('Global C-alpha RMSD (Å)')
ax2.set_xlabel('Loop Length (#Residues)')
ax2.plot(x, b,'o',markersize=1)
ax2.text(0,20, "Pearson correlation coefficient:{}".format(round(r2[0],2)), fontsize=10, fontweight='bold')
ax2.plot(np.unique(x), m2*np.unique(x) + b2, label="y={}x+{}".format(round(m2, 2), round(b2, 2)))
ax2.legend(loc="center left", bbox_to_anchor=(1, 0.5))

plt.show()
fig.savefig("RMSDvsLoopLength_graph.png")