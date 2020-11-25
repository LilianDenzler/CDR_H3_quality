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
#GRAPH: CDRH3 Loop Ca-RMSD vs Loop Length -violin plot !!
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

#separate values by loop length
dic1={}
values1=[]
averages=[]
list1=np.arange(1, max(x), 4)
for q, a in zip(x, a):
	for i in range (len(list1)-1):
		if int(q)>list1[i] and int(q)<=list1[i+1]:
			average=(list1[i]+list1[i+1])/2
			if average in dic1.keys():
				dic1[average].append(a)
			else:
				dic1.update({average:[a]})
				averages.append(average)
[values1.append(a) for q, a in dic1.items() ]

q1=[]
med=[]
q3=[]
for i in values1:
	quartile1, medians, quartile3 = np.percentile(i,[25, 50, 75])
	q1.append(quartile1)
	q3.append(quartile3)
	med.append(medians)

dic2={}
values2=[]
list2=np.arange(1, max(x), 4)
for q, b in zip(x, b):
	for i in range (len(list2)-1):
		if int(q)>list2[i] and int(q)<=list2[i+1]:
			average=(list2[i]+list2[i+1])/2
			if average in dic2.keys():
				dic2[average].append(b)
			else:
				dic2.update({average:[b]})
[values2.append(b) for q, b in dic2.items() ]

q1_2=[]
med_2=[]
q3_2=[]
for i in values2:
	quartile1_2, medians_2, quartile3_2 = np.percentile(i,[25, 50, 75])
	q1_2.append(quartile1_2)
	q3_2.append(quartile3_2)
	med_2.append(medians_2)

fig, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(50, 50))# sharex=True, sharey=True
plt.subplots_adjust(wspace=0.5, hspace=0.5)
fig.suptitle('Loop Length vs CDRH3 Loop RMSD', size="20", fontweight='bold')

ax1.set_ylabel('Local C-alpha RMSD (Å)')
ax1.set_xlabel('Loop Length (#Residues)')
#quartile1, medians, quartile3 = np.percentile(values1,[25, 50, 75])
ax1.violinplot(values1, dic1.keys(),  widths=2, showmeans=True, showmedians=False)
ax1.vlines(dic1.keys(), q1, q3, color='k', linestyle='-', lw=5)
ax1.text(1,6.5, "Pearson correlation coefficient:{}".format(round(r1[0],2)), fontsize=10, fontweight='bold')
ax1.set_xticks(list1)
#ax1.plot(np.unique(x), m1*np.unique(x) + b1, label="y={}x+{}".format(round(m1, 2), round(b1, 2)))
ax1.scatter(dic1.keys(), med, color='white', s=5, zorder=3)
ax1.legend(loc="center left", bbox_to_anchor=(1, 0.5))

ax2.set_ylabel('Global C-alpha RMSD (Å)')
ax2.set_xlabel('Loop Length (#Residues)')
ax2.violinplot(values2, dic2.keys(),  widths=1.6, showmeans=True, showmedians=False)
ax2.vlines(dic2.keys(), q1_2, q3_2, color='k', linestyle='-', lw=5)
ax2.scatter(dic2.keys(), med_2, color='white', s=5, zorder=3)
ax2.text(1,20, "Pearson correlation coefficient:{}".format(round(r2[0],2)), fontsize=10, fontweight='bold')
ax2.set_xticks(list2)

#ax2.plot(np.unique(x), m2*np.unique(x) + b2, label="y={}x+{}".format(round(m2, 2), round(b2, 2)))
ax2.legend(loc="center left", bbox_to_anchor=(1, 0.5))
#ax2.set_title("local_CA vs global_CA", fontsize=14, fontweight='bold')


plt.show()
fig.savefig("length_vs_RMSD_violin.png")
