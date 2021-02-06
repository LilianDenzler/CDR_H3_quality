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
def run(threshold):
	data = pandas.read_csv(sys.argv[1])
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
	print(data.shape)
	for i in range(data.shape[0]):
		#print(i)
		row=data.iloc[i]
		c+=1
		#print(row)
		if row[0]=="local_AA":
			pass
		if row[0]=="local_CA":
			pass
		if row[0]=="global_AA":
			pass
		if row[0]=="global_CA":
			pass
		else:
			key = str(data.iloc[i]['length'])
			#print(key)
			if key not in result_lAA:
				#print(key)
				result_lAA.update({key: [data.iloc[i]['local_AA']]})
				result_lCA.update({key: [data.iloc[i]['local_CA']]})
				result_gAA.update({key: [data.iloc[i]['global_AA']]})
				result_gCA.update({key: [data.iloc[i]['global_CA']]})
			else:
				result_lAA[key].append(str(data.iloc[i]['local_AA']))
				result_lCA[key].append(str(data.iloc[i]['local_CA']))
				result_gAA[key].append(str(data.iloc[i]['global_AA']))
				result_gCA[key].append(str(data.iloc[i]['global_CA']))

	print(len(result_lAA.keys()), c)
	#print(result_lAA)

	for asd in [result_lAA,result_gCA,result_lCA, result_gAA]:
		for keys,values in asd.items():
			#print(keys,values)
			asd[keys]=[float(i) for i in asd[keys]]
			for n,i in enumerate(asd[keys]):
				if i< threshold:
					#print("hey")
					asd[keys][n]=1
				else:
					asd[keys][n]=0
			#print("keys", a[keys])
			asd[keys]=sum(asd[keys])/len(asd[keys])
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
	print(r1, r2, r3, r4)


	return(x_lAA,y_lAA,x_gAA, y_gAA, x_lCA, y_lCA, x_gCA, y_gCA)

x_lAA,y_lAA,x_gAA, y_gAA, x_lCA, y_lCA, x_gCA, y_gCA= run(2)
ax1=sns.scatterplot(x=x_lAA, y=y_lAA, s=50, color="#FF00FF", label="< 2Å")
x_lAA,y_lAA,x_gAA, y_gAA, x_lCA, y_lCA, x_gCA, y_gCA= run(3)
ax1=sns.scatterplot(x=x_lAA, y=y_lAA, s=50, color="#8200FF", label="< 3Å")
x_lAA,y_lAA,x_gAA, y_gAA, x_lCA, y_lCA, x_gCA, y_gCA= run(4)
ax1=sns.scatterplot(x=x_lAA, y=y_lAA, s=50, color="#0000BC", label="< 4Å")
ax1.set_xlabel('Loop Length l')
ax1.set_ylabel('% of Structures (with length=l) below threshold')
ax1.set_title("Performance for different loop lengths (Local atom-atom)", fontsize=14, fontweight='bold')
ax1.legend(loc="center right")
#plt.xlim([-0.6,1.01])
#plt.ylim([-0.1,8.5])
plt.savefig("local_AA_length_smaller2.png")
plt.show()


x_lAA,y_lAA,x_gAA, y_gAA, x_lCA, y_lCA, x_gCA, y_gCA= run(2)
ax1=sns.scatterplot(x=x_lCA, y=y_lCA, s=50, color="#FF00FF", label="< 2Å")
x_lAA,y_lAA,x_gAA, y_gAA, x_lCA, y_lCA, x_gCA, y_gCA= run(3)
ax1=sns.scatterplot(x=x_lCA, y=y_lCA, s=50, color="#8200FF", label="< 3Å")
x_lAA,y_lAA,x_gAA, y_gAA, x_lCA, y_lCA, x_gCA, y_gCA= run(4)
ax1=sns.scatterplot(x=x_lCA, y=y_lCA, s=50, color="#0000BC", label="< 4Å")
ax1.set_xlabel('Loop Length l')
ax1.set_ylabel('% of Structures (with length=l) below threshold')
ax1.legend(loc="center right")
#plt.xlim([-0.6,1.01])
#plt.ylim([-0.1,8.5])
ax1.set_title("Performance for different loop lengths (Local C-alpha)", fontsize=14, fontweight='bold')
plt.savefig("local_CA_length_smaller2.png")
plt.show()

x_lAA,y_lAA,x_gAA, y_gAA, x_lCA, y_lCA, x_gCA, y_gCA= run(2)
ax1=sns.scatterplot(x=x_gAA, y=y_gAA, s=50, color="#FF00FF", label="< 2Å")
x_lAA,y_lAA,x_gAA, y_gAA, x_lCA, y_lCA, x_gCA, y_gCA= run(3)
ax1=sns.scatterplot(x=x_gAA, y=y_gAA, s=50, color="#8200FF", label="< 3Å")
x_lAA,y_lAA,x_gAA, y_gAA, x_lCA, y_lCA, x_gCA, y_gCA= run(4)
ax1=sns.scatterplot(x=x_gAA, y=y_gAA, s=50, color="#0000BC", label="< 4Å")
ax1.set_xlabel('Loop Length l')
ax1.set_ylabel('% of Structures (with length=l) below threshold')
ax1.set_title("Performance for different loop lengths (Global atom-atom)", fontsize=14, fontweight='bold')
ax1.legend(loc="center right")
#plt.xlim([-0.6,1.01])
#plt.ylim([-0.1,8.5])
plt.savefig("global_AA_length_smaller2.png")
plt.show()


x_lAA,y_lAA,x_gAA, y_gAA, x_lCA, y_lCA, x_gCA, y_gCA= run(2)
ax1=sns.scatterplot(x=x_gCA, y=y_gCA, s=50, color="#FF00FF", label="< 2Å")
x_lAA,y_lAA,x_gAA, y_gAA, x_lCA, y_lCA, x_gCA, y_gCA= run(3)
ax1=sns.scatterplot(x=x_gCA, y=y_gCA, s=50, color="#8200FF", label="< 3Å")
x_lAA,y_lAA,x_gAA, y_gAA, x_lCA, y_lCA, x_gCA, y_gCA= run(4)
ax1=sns.scatterplot(x=x_gCA, y=y_gCA, s=50, color="#0000BC", label="< 4Å")
ax1.set_xlabel('Loop Length l')
ax1.set_ylabel('% of Structures (with length=l) below threshold')
ax1.set_title("Performance for different loop lengths (Global C-alpha)", fontsize=14, fontweight='bold')
ax1.legend(loc="center right")
#plt.xlim([-0.6,1.01])
#plt.ylim([-0.1,8.5])
plt.savefig("global_CA_length_smaller2.png")
plt.show()

'''
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=False, sharey=False, figsize=(50, 50))# sharex=True, sharey=True
plt.subplots_adjust(wspace=0.5, hspace=0.5)
fig.suptitle('Loop length vs CDRH3 RMSD', size="20", fontweight='bold')
ax1.set_xlabel('Loop Length l')
ax1.set_ylabel('% of Structures (with length=l) with < 2Å')
ax1.plot(x_lAA, y_lAA,'o')
#ax1.plot(x_lAA, new_y, label="fit");
#ax1.legend(loc='best')
ax1.text(0,15, "Pearson Correlation coefficient:{}".format(round(r1[0],2)), fontsize=10, fontweight='bold')
ax1.set_title("Local Atom", fontsize=14, fontweight='bold')

ax2.set_xlabel('loop length l')
ax2.set_ylabel('% of structures (with length=l) with < 2Å')
ax2.plot(x_gAA, y_gAA,'o',markersize=5)
ax2.text(0,15, "Pearson Correlation coefficient:{}".format(round(r2[0],2)), fontsize=10, fontweight='bold')
#ax2.plot(np.unique(x_gAA), m2*np.unique(x_gAA) + b2, label="y={}x+{}".format(round(m2, 2), round(b2, 2)))
#ax2.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax2.set_title("Global ATom", fontsize=14, fontweight='bold')

ax3.set_xlabel('loop length l')
ax3.set_ylabel('% of structures (with length=l) with < 2Å')
ax3.plot(x_lCA, y_lCA, 'o',markersize=5)
ax3.text(0,5, "Pearson Correlation coefficient:{}".format(round(r3[0],2)), fontsize=10, fontweight='bold')
#ax3.plot(np.unique(x_lCA), m3*np.unique(x_lCA) + b3, label="y={}x+{}".format(round(m3, 2), round(b3, 2)))
#ax3.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax3.set_title("Local C-alpha", fontsize=14, fontweight='bold')

ax4.set_xlabel('loop length l')
ax4.set_ylabel('% of structures (with length=l) with < 2Å')
ax4.plot(x_gCA, y_gCA, 'o',markersize=5)
ax4.text(0,15, "Pearson Correlation coefficient:{}".format(round(r4[0],2)), fontsize=10, fontweight='bold')
#ax4.plot(np.unique(x_gCA), m4*np.unique(x_gCA) + b4, label="y={}x+{}".format(round(m4, 2), round(b4, 2)))
#ax4.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax4.set_title("Global C-alpha", fontsize=14, fontweight='bold')

plt.show()
fig.savefig("RMSD_graph_2Angstrom.png")
'''