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
#GRAPH: Ca-RMSD under 2A vs Loop Length
#DO FITTING!!!!!!!!!!!!!!!!
#################################################################
columns=["local_AA","local_CA","global_AA","global_CA","ID", "length", "identity", "similarity", "template", "target"]
data = pandas.read_csv(sys.argv[1], names=columns)
x=data.local_AA.tolist()
y=data.global_AA.tolist()
a=data.local_CA.tolist()
b=data.global_CA.tolist()
c=data.length.tolist()
del x[0]
del y[0]
del a[0]
del b[0]
del c[0]
total=len(x)
x = [g for g in x if str(g) != 'nan']
x = [float(g) for g in x]
y = [g for g in y if str(g) != 'nan']
y = [float(g) for g in y]
a = [g for g in a if str(g) != 'nan']
a = [float(g) for g in a]
b = [g for g in b if str(g) != 'nan']
b = [float(g) for g in b]
c= [g for g in c if str(g) !='nan']
c= [float(g) for g in c]

#get number of under 2A 
reader = csv.reader(open(sys.argv[1]))

result_lAA = {}
result_gAA = {}
result_lCA = {}
result_gCA = {}
counter_lAA=0
counter_gAA=0
counter_lCA=0
counter_gCA=0
c=0
for row in reader:
	c+=1
	#print(row)
	if row[0]=="local_AA":
		pass
	else:
		key = row[5]
		if key not in result_lAA:
			result_lAA.update({key: [row[0]]})
			result_lCA.update({key: [row[1]]})
			result_gAA.update({key: [row[2]]})
			result_gCA.update({key: [row[3]]})
		else:
			result_lAA[key].append(row[0])
			result_lCA[key].append(row[1])
			result_gAA[key].append(row[2])
			result_gCA[key].append(row[3])

print(len(result_lAA.keys()), c)
#print(result_lAA)

for a in [result_lAA,result_gCA,result_lCA, result_gAA]:
	for keys,values in a.items():
		a[keys]=[float(i) for i in a[keys]]
		for n,i in enumerate(a[keys]):
			if i< 2:
				#print("hey")
				a[keys][n]=1
			else:
				a[keys][n]=0
		#print("keys", a[keys])
		a[keys]=sum(a[keys])/len(a[keys])
		#a = {float(k):float(v) for k,v in a.items()}



lists = (result_lAA.items()) 
x_lAA, y_lAA = zip(*lists) 
x_lAA=list(x_lAA)
y_lAA=list(y_lAA)
x_lAA = [int(v) for v in x_lAA]

lists = (result_gAA.items()) 
x_gAA, y_gAA = zip(*lists) 
x_gAA=list(x_gAA)
y_gAA=list(y_gAA)
x_gAA = [int(v) for v in x_gAA]

lists = (result_lCA.items())
x_lCA, y_lCA = zip(*lists)
x_lCA=list(x_lCA)
y_lCA=list(y_lCA)
x_lCA = [int(v) for v in x_lCA] 

lists = (result_gCA.items()) 
x_gCA, y_gCA = zip(*lists)
x_gCA=list(x_gCA)
y_gCA=list(y_gCA)
x_gCA = [int(v) for v in x_gCA]


r1=stats.pearsonr(x_lAA, y_lAA) #	(Pearson’s correlation coefficient, 2-tailed p-value)
r2=stats.pearsonr(x_gAA, y_gAA)
r3=stats.pearsonr(x_lCA, y_lCA)
r4=stats.pearsonr(x_gCA, y_gCA)
x_lAA=np.array(x_lAA)
y_lAA=np.array(y_lAA)
print(x_lAA)
print(type(y_lAA[1]))

def test_func(x, a, b, c):
     return a * np.exp(-b * x) + c

params, params_covariance = optimize.curve_fit(test_func, x_lAA, y_lAA)
a, b, c= params
new_y=test_func(x_lAA, a,b,c)
print(r1, r2, r3, r4)
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=False, sharey=False, figsize=(50, 50))# sharex=True, sharey=True
plt.subplots_adjust(wspace=0.5, hspace=0.5)
fig.suptitle('Loop length vs CDRH3 RMSD', size="20", fontweight='bold')
ax1.set_xlabel('loop length l')
ax1.set_ylabel('% of structures (with length=l) with > 2Å')
ax1.plot(x_lAA, y_lAA,'o')
#ax1.plot(x_lAA, new_y, label="fit");
#ax1.legend(loc='best')
ax1.text(0,15, "Pearson Correlation coefficient:{}".format(round(r1[0],2)), fontsize=10, fontweight='bold')
ax1.set_title("Local Atom", fontsize=14, fontweight='bold')

ax2.set_xlabel('loop length l')
ax2.set_ylabel('% of structures (with length=l) with > 2Å')
ax2.plot(x_gAA, y_gAA,'o',markersize=5)
ax2.text(0,15, "Pearson Correlation coefficient:{}".format(round(r2[0],2)), fontsize=10, fontweight='bold')
#ax2.plot(np.unique(x_gAA), m2*np.unique(x_gAA) + b2, label="y={}x+{}".format(round(m2, 2), round(b2, 2)))
#ax2.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax2.set_title("Global ATom", fontsize=14, fontweight='bold')

ax3.set_xlabel('loop length l')
ax3.set_ylabel('% of structures (with length=l) with > 2Å')
ax3.plot(x_lCA, y_lCA, 'o',markersize=5)
ax3.text(0,5, "Pearson Correlation coefficient:{}".format(round(r3[0],2)), fontsize=10, fontweight='bold')
#ax3.plot(np.unique(x_lCA), m3*np.unique(x_lCA) + b3, label="y={}x+{}".format(round(m3, 2), round(b3, 2)))
#ax3.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax3.set_title("Local C-alpha", fontsize=14, fontweight='bold')

ax4.set_xlabel('loop length l')
ax4.set_ylabel('% of structures (with length=l) with > 2Å')
ax4.plot(x_gCA, y_gCA, 'o',markersize=5)
ax4.text(0,15, "Pearson Correlation coefficient:{}".format(round(r4[0],2)), fontsize=10, fontweight='bold')
#ax4.plot(np.unique(x_gCA), m4*np.unique(x_gCA) + b4, label="y={}x+{}".format(round(m4, 2), round(b4, 2)))
#ax4.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax4.set_title("Global C-alpha", fontsize=14, fontweight='bold')

plt.show()
fig.savefig("RMSD_graph_2Angstrom.png")
