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
#GRAPH: Sequence Similarity/Identity vs RMSD
#DO FITTING!!!!!!!!!!!!!!!!
#################################################################
columns=["local_AA","local_CA","global_AA","global_CA","ID", "length", "identity", "similarity", "template", "target"]
data = pandas.read_csv(sys.argv[1], names=columns)
x=data.identity.tolist()
y=data.similarity.tolist()
a=data.local_CA.tolist()

del x[0]
del y[0]
del a[0]

x = [g for g in x if str(g) != 'nan']
x = [float(g) for g in x]
y = [g for g in y if str(g) != 'nan']
y = [float(g) for g in y]
a = [g for g in a if str(g) != 'nan']
a = [float(g) for g in a]
print(len(x),len(a))
r1=stats.pearsonr(x, a) #	(Pearson’s correlation coefficient, 2-tailed p-value)
r2=stats.pearsonr(y,a)

def func(x, a, b):
	return a*np.power(x, b)

m1, b1 = np.polyfit(x, a, 1)
m2, b2 = np.polyfit(y, a, 1)

print(r1, r2)
#fig, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(50, 50))# sharex=True, sharey=True#
pl=sns.jointplot(x,a)
pl.title("Sequence Identity ", fontsize=14, fontweight='bold')
pl.xlabel('seq identity') 
pl.ylabel('RMSD Local CA (Å)')

pl.text(0,7, "Pearson Correlation coefficient:{}".format(round(r1[0],2)), fontsize=10, fontweight='bold')


ax2.set_xlabel('seq similarity')
ax2.set_ylabel('RMSD Local CA (Å)')
ax2.plot(y, a,'o',markersize=5)
ax2.text(0,8, "Pearson Correlation coefficient:{}".format(round(r2[0],2)), fontsize=10, fontweight='bold')
ax2.set_title("Sequence Similarity", fontsize=14, fontweight='bold')

plt.show()
fig.savefig("seqId_similarity_RMSD.png") 