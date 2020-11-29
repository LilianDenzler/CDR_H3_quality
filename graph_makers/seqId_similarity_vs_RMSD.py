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
from itertools import chain

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
print(len(x),len(a))
r1=stats.pearsonr(x, a) #	(Pearson’s correlation coefficient, 2-tailed p-value)
r2=stats.pearsonr(y,a)

print(r1, r2)
df1 = pd.DataFrame(list(zip(x, a)), columns =['sequence identity', 'RMSD Local CA (Å)']) 
df2 = pd.DataFrame(list(zip(y, a)), columns =['sequence similarity', 'RMSD Local CA (Å)']) 
sns.jointplot(data=df1, x="sequence identity", y="RMSD Local CA (Å)", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.xlim([-0.05,1.01])
plt.ylim([-0.1,8.5])
sns.jointplot(data=df2, x="sequence similarity", y="RMSD Local CA (Å)", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.xlim([-0.6,1.01])
plt.ylim([-0.1,8.5])

plt.show()
#plt.savefig("seqId_similarity_RMSD.png") 