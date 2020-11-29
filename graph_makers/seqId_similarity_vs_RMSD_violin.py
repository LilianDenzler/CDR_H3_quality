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
columns=["local_AA","local_CA","global_AA","global_CA","ID", "length", "identity", "similarity", "template", "target"]
data = pd.read_csv(sys.argv[1], names=columns)
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

r1=stats.pearsonr(x, a) #	(Pearson’s correlation coefficient, 2-tailed p-value)
r2=stats.pearsonr(y,a)

list1=[x,a]
df1 = pd.DataFrame(list1)
df1=pd.DataFrame.transpose(df1)
df1.rename(columns={0: 'identity', 1: 'localCA'}, inplace=True)
print(df1)
fig, axes = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(50, 50))# sharex=True, sharey=True
axes[0].set_xlabel('seq identity')
axes[0].set_ylabel('RMSD LocalCA (Å)')
axes[0].plot(x, a,'o',markersize=5)
axes[0].text(0,22, "Pearson Correlation coefficient:{}".format(round(r1[0],2)), fontsize=10, fontweight='bold')
axes[0].set_title("Sequence Identity ", fontsize=14, fontweight='bold')

axes[1].set_xlabel('seq similarity')
axes[1].set_ylabel('RMSD LocalCA (Å)')
axes[1].plot(y, a,'o',markersize=5)
axes[1].text(-0.4,22, "Pearson Correlation coefficient:{}".format(round(r2[0],2)), fontsize=10, fontweight='bold')
axes[1].set_title("Sequence Similarity", fontsize=14, fontweight='bold')

plt.show()
fig.savefig("seqId_similarity_RMSD.png") 