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

r1=stats.pearsonr(x, y) #	(Pearson’s correlation coefficient, 2-tailed p-value)

fig = plt.figure()

plt.xlabel('Hydrophobicity')
plt.ylabel('RMSD Local CA (Å)')
plt.plot(x, y,'o',markersize=5)
plt.text(-10,8, "Pearson Correlation coefficient:{}".format(round(r1[0],2)), fontsize=10, fontweight='bold')
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
sns.jointplot(x=x,y=y,kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.show()
fig.savefig("Hydrophobicity_diff_RMSD.png") 
