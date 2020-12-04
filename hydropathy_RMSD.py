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
#GRAPH: Ca-RMSD vs Hydrophobicity
#DO FITTING!!!!!!!!!!!!!!!!
#################################################################
data = pd.read_csv(sys.argv[1], header=0)
x=data.Hydropathy.tolist()
y=data.local_CA.tolist()


del x[0]
del y[0]
x = [float(k) for k in x]
y = [float(k) for k in y]
print(type(x[1]))
print(type(y[1]))

r1=stats.pearsonr(x, y) #	(Pearson’s correlation coefficient, 2-tailed p-value)

print(r1)
fig = plt.figure()
plt.xlabel('Hydrophobicity')
plt.ylabel('RMSD Local CA (Å)')
plt.plot(x, y,'o',markersize=5)
plt.text(-10,8, "Pearson Correlation coefficient:{}".format(round(r1[0],2)), fontsize=10, fontweight='bold')
#plt.plot(np.unique(x), m1*np.unique(x) + b1, label="y={}x+{}".format(round(m1, 2), round(b1, 2)))
#popt, pcov = curve_fit(func, np.asarray(x), np.asarray(y))
#plt.plot(x, func(np.asarray(x), *popt),label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
#plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.title("Hydrophobicity vs RMSD", fontsize=14, fontweight='bold')
plt.ylim([0,9])
plt.show()
fig.savefig("Hydrophobicity_RMSD.png") 

x=data.Hydropathy_diff.tolist()
y=data.local_CA.tolist()


del x[0]
del y[0]
x = [float(k) for k in x]
y = [float(k) for k in y]
print(type(x[1]))
print(type(y[1]))

r1=stats.pearsonr(x, y) #	(Pearson’s correlation coefficient, 2-tailed p-value)


print(r1)
fig = plt.figure()
sns.jointplot(x,y,kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
#plt.text(-40,22, "Pearson Correlation coefficient:{}".format(round(r1[0],2)), fontsize=10, fontweight='bold')
#plt.plot(np.unique(x), m1*np.unique(x) + b1, label="y={}x+{}".format(round(m1, 2), round(b1, 2)))
#popt, pcov = curve_fit(func, np.asarray(x), np.asarray(y))
#plt.plot(x, func(np.asarray(x), *popt),label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
#plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
#plt.title("Hydrophobicity difference to template vs RMSD", fontsize=14, fontweight='bold')
plt.ylim([0,9])
plt.xlim([-0.05,2])
plt.show()
fig.savefig("Hydrophobicity_diff_RMSD.png") 



