#!/usr/bin/env python
#script for graphing the different RMSD values against each other
import sys
import os
import pandas
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
columns=["local_AA","local_CA","global_AA","global_CA","ID"]
data = pandas.read_csv(sys.argv[1], names=columns)
x=data.local_AA.tolist()
y=data.global_AA.tolist()
a=data.local_CA.tolist()
b=data.global_CA.tolist()
del x[0]
del y[0]
del a[0]
del b[0]
x = [g for g in x if str(g) != 'nan']
x = [float(g) for g in x]
y = [g for g in y if str(g) != 'nan']
y = [float(g) for g in y]
a = [g for g in a if str(g) != 'nan']
a = [float(g) for g in a]
b = [g for g in b if str(g) != 'nan']
b = [float(g) for g in b]

r1=stats.pearsonr(x, y) #	(Pearson’s correlation coefficient, 2-tailed p-value)
r2=stats.pearsonr(a, b)
r3=stats.pearsonr(x, a)
r4=stats.pearsonr(y, b)

m1, b1 = np.polyfit(x, y, 1)
m2, b2 = np.polyfit(a, b, 1)
m3, b3 = np.polyfit(x, a, 1)
m4, b4 = np.polyfit(y, b, 1)
print(r1, r2, r3, r4)
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=False, sharey=False, figsize=(50, 50))# sharex=True, sharey=True
plt.subplots_adjust(wspace=0.5, hspace=0.5)
fig.suptitle('Comparison of CDRH3 loop RMSD Calculations', size="20", fontweight='bold')
ax1.set_xlabel('local_AA')
ax1.set_ylabel('global_AA')
ax1.plot(x, y,'o',markersize=1)
ax1.text(0,15, "Pearson Correlation coefficient:{}".format(round(r1[0],2)), fontsize=10, fontweight='bold')
ax1.plot(np.unique(x), m1*np.unique(x) + b1, label="y={}x+{}".format(round(m1, 2), round(b1, 2)))
ax1.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax1.set_title("local_AA vs global_AA", fontsize=14, fontweight='bold')
ax2.set_xlabel('local_CA')
ax2.set_ylabel('global_CA')
ax2.plot(a, b,'o',markersize=1)
ax2.text(0,15, "Pearson Correlation coefficient:{}".format(round(r2[0],2)), fontsize=10, fontweight='bold')
ax2.plot(np.unique(a), m2*np.unique(a) + b2, label="y={}x+{}".format(round(m2, 2), round(b2, 2)))
ax2.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax2.set_title("local_CA vs global_CA", fontsize=14, fontweight='bold')
ax3.set_xlabel('local_AA')
ax3.set_ylabel('local_CA')
ax3.plot(x, a, 'o',markersize=1)
ax3.text(0,5, "Pearson Correlation coefficient:{}".format(round(r3[0],2)), fontsize=10, fontweight='bold')
ax3.plot(np.unique(x), m3*np.unique(x) + b3, label="y={}x+{}".format(round(m3, 2), round(b3, 2)))
ax3.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax3.set_title("local_AA vs local_CA", fontsize=14, fontweight='bold')
ax4.set_xlabel('global_AA')
ax4.set_ylabel('global_CA')
ax4.plot(y, b, 'o',markersize=1)
ax4.text(0,15, "Pearson Correlation coefficient:{}".format(round(r4[0],2)), fontsize=10, fontweight='bold')
ax4.plot(np.unique(y), m4*np.unique(y) + b4, label="y={}x+{}".format(round(m4, 2), round(b4, 2)))
ax4.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax4.set_title("global_AA vs global_CA", fontsize=14, fontweight='bold')

plt.show()
fig.savefig("RMSD_graph.png")
