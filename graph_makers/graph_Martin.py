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
'''
#################################################################
#GRAPH: Ca-RMSD under 2A vs seqId and similarity
#DO FITTING!!!!!!!!!!!!!!!!
#################################################################
columns=["local_AA","local_CA","global_AA","global_CA","ID", "length", "identity", "similarity", "template", "target"]
data = pandas.read_csv(sys.argv[1], names=columns)
a=data.global_AA.tolist()
x=data.similarity.tolist()
y=data.identity.tolist()
del x[0]
del y[0]
del a[0]
total=len(a)


#get number of under 2A 
reader = csv.reader(open(sys.argv[1]))
###################
############
###########
##############
###########

def func(x, a, b):
	return a*np.power(x, b)

#print(r1, r2)
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, sharey=False, figsize=(50, 50))# sharex=True, sharey=True
plt.subplots_adjust(wspace=0.5, hspace=0.5)
#fig.suptitle('CDRH3 loop RMSD and Loop Length', size="20", fontweight='bold')
ax1.set_xlabel('Similarity')
ax1.set_ylabel('')
ax1.plot(x_gAA, y_gAA,'o',markersize=5)
ax1.plot(xb_gAA, yb_gAA,'ro',markersize=5)
#ax1.text(0,15, "Pearson Correlation coefficient:{}".format(round(r1[0],2)), fontsize=10, fontweight='bold')
#ax1.plot(np.unique(x_lAA), m1*np.unique(x_lAA) + b1, label="y={}x+{}".format(round(m1, 2), round(b1, 2)))
#popt, pcov = curve_fit(func, np.asarray(x_lAA), np.asarray(y_lAA))
#ax1.plot(x_lAA, func(np.asarray(x_lAA), *popt),label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
#ax1.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax1.set_title("Similarity vs RMSD under 2Å", fontsize=14, fontweight='bold')

ax2.set_xlabel('Identity')
ax2.set_ylabel('% of structures with > 2Å')
ax2.plot(x2_gAA, y2_gAA,'o',markersize=5)
#ax2.text(0,15, "Pearson Correlation coefficient:{}".format(round(r2[0],2)), fontsize=10, fontweight='bold')
#ax2.plot(np.unique(x_gAA), m2*np.unique(x_gAA) + b2, label="y={}x+{}".format(round(m2, 2), round(b2, 2)))
#ax2.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax2.set_title("Identity vs RMSD under 2Å", fontsize=14, fontweight='bold')


plt.show()
fig.savefig("seqId_similarity_RMSD_2A.png")
#################################################################
#GRAPH: Ca-RMSD vs Hydrophobicity
#DO FITTING!!!!!!!!!!!!!!!!
#################################################################
columns=["charge", "ID", "Surface", "Hydropathy", "local_AA","local_CA","global_AA","global_CA"]
data = pandas.read_csv(sys.argv[1], names=columns)
x=data.Hydropathy.tolist()
y=data.global_AA.tolist()


del x[0]
del y[0]
x = [float(k) for k in x]
y = [float(k) for k in y]
print(type(x[1]))
print(type(y[1]))

r1=stats.pearsonr(x, y) #	(Pearson’s correlation coefficient, 2-tailed p-value)


def func(x, a, b):
	return a*np.power(x, b)

m1, b1 = np.polyfit(x, y, 1)

print(r1)
fig = plt.figure()
plt.xlabel('Hydrophobicity')
plt.ylabel('RMSD Global AA (Å)')
plt.plot(x, y,'o',markersize=5)
plt.text(-40,22, "Pearson Correlation coefficient:{}".format(round(r1[0],2)), fontsize=10, fontweight='bold')
#plt.plot(np.unique(x), m1*np.unique(x) + b1, label="y={}x+{}".format(round(m1, 2), round(b1, 2)))
popt, pcov = curve_fit(func, np.asarray(x), np.asarray(y))
plt.plot(x, func(np.asarray(x), *popt),label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
#plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.title("Hydropathy", fontsize=14, fontweight='bold')
plt.show()
fig.savefig("Hydropathy_RMSD.png") 


#################################################################
#GRAPH: Ca-RMSD vs Surface
#DO FITTING!!!!!!!!!!!!!!!!
#################################################################
columns=["charge", "ID", "Surface", "Hydropathy", "local_AA","local_CA","global_AA","global_CA"]
data = pandas.read_csv(sys.argv[1], names=columns)
x=data.Surface.tolist()
y=data.global_AA.tolist()


del x[0]
del y[0]
x = [float(k) for k in x]
y = [float(k) for k in y]


bady=[]
badx=[]
goodx=[]
goody=[]
for i in range (len(y)):
	if y[i] > 3:
		bady.append(y[i])
		badx.append(x[i])
	else:
		goody.append(y[i])
		goodx.append(x[i])
r1=stats.pearsonr(goodx, goody) #	(Pearson’s correlation coefficient, 2-tailed p-value)
r2=stats.pearsonr(badx, bady)

m1, b1 = np.polyfit(goodx, goody, 1)
m2, b2 = np.polyfit(badx, bady, 1)

print(r1)
fig = plt.figure()
plt.xlabel('Surface Accessibility')
plt.ylabel('RMSD Global AA (Å)')
plt.plot(goodx, goody,'o',markersize=5)
plt.plot(badx, bady,'ro',markersize=5)
plt.text(-40,20, "Pearson Correlation coefficient:{}".format(round(r1[0],2)), fontsize=10, fontweight='bold', color="blue")
plt.text(-40,22, "Pearson Correlation coefficient:{}".format(round(r2[0],2)), fontsize=10, fontweight='bold', color="red")
plt.plot(np.unique(goodx), m1*np.unique(goodx) + b1, label="y={}x+{}".format(round(m1, 2), round(b1, 2)))
plt.plot(np.unique(badx), m2*np.unique(badx) + b2, label="y={}x+{}".format(round(m2, 2), round(b2, 2)))
#popt, pcov = curve_fit(func, np.asarray(x), np.asarray(y))
#plt.plot(x, func(np.asarray(x), *popt),label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
#plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.title("Surface Accessibility", fontsize=14, fontweight='bold')
plt.show()
fig.savefig("surface_RMSD.png") 
###########################

columns=["charge", "ID", "Surface", "Hydropathy", "local_AA","local_CA","global_AA","global_CA"]
data = pandas.read_csv(sys.argv[1], names=columns)
x=data.Surface.tolist()
y=data.local_CA.tolist()


del x[0]
del y[0]
x = [float(k) for k in x]
y = [float(k) for k in y]


bady=[]
badx=[]
goodx=[]
goody=[]
for i in range (len(y)):
	if y[i] > 2:
		bady.append(y[i])
		badx.append(x[i])
	else:
		goody.append(y[i])
		goodx.append(x[i])
r1=stats.pearsonr(goodx, goody) #	(Pearson’s correlation coefficient, 2-tailed p-value)
r2=stats.pearsonr(badx, bady)

m1, b1 = np.polyfit(goodx, goody, 1)
m2, b2 = np.polyfit(badx, bady, 1)

print(r1)
fig = plt.figure()
plt.xlabel('Surface Accessibility')
plt.ylabel('RMSD Local CA (Å)')
plt.plot(goodx, goody,'o',markersize=5)
plt.plot(badx, bady,'ro',markersize=5)
plt.text(-40,7, "Pearson Correlation coefficient:{}".format(round(r1[0],2)), fontsize=10, fontweight='bold', color="blue")
plt.text(-40,8, "Pearson Correlation coefficient:{}".format(round(r2[0],2)), fontsize=10, fontweight='bold', color="red")
plt.plot(np.unique(goodx), m1*np.unique(goodx) + b1, label="y={}x+{}".format(round(m1, 2), round(b1, 2)))
plt.plot(np.unique(badx), m2*np.unique(badx) + b2, label="y={}x+{}".format(round(m2, 2), round(b2, 2)))
#popt, pcov = curve_fit(func, np.asarray(x), np.asarray(y))
#plt.plot(x, func(np.asarray(x), *popt),label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
#plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.title("Surface Accessibility", fontsize=14, fontweight='bold')
plt.show()
fig.savefig("surface_RMSD_lCA.png") 

#################################################################
#GRAPH: Ca-RMSD vs Surface
#DO FITTING!!!!!!!!!!!!!!!!
#################################################################
columns=["charge", "ID", "Surface", "Hydropathy", "local_AA","local_CA","global_AA","global_CA"]
data = pandas.read_csv(sys.argv[1], names=columns)
x=data.charge.tolist()
y=data.global_AA.tolist()


del x[0]
del y[0]
x = [float(k) for k in x]
y = [float(k) for k in y]


bady=[]
badx=[]
goodx=[]
goody=[]
for i in range (len(y)):
	if y[i] > 3:
		bady.append(y[i])
		badx.append(x[i])
	else:
		goody.append(y[i])
		goodx.append(x[i])
r1=stats.pearsonr(goodx, goody) #	(Pearson’s correlation coefficient, 2-tailed p-value)
r2=stats.pearsonr(badx, bady)

m1, b1 = np.polyfit(goodx, goody, 1)
m2, b2 = np.polyfit(badx, bady, 1)

print(r1)
fig = plt.figure()
plt.xlabel('Surface Accessibility')
plt.ylabel('RMSD Global AA (Å)')
plt.plot(goodx, goody,'o',markersize=5)
plt.plot(badx, bady,'ro',markersize=5)
plt.text(-40,20, "Pearson Correlation coefficient:{}".format(round(r1[0],2)), fontsize=10, fontweight='bold', color="blue")
plt.text(-40,22, "Pearson Correlation coefficient:{}".format(round(r2[0],2)), fontsize=10, fontweight='bold', color="red")
plt.plot(np.unique(goodx), m1*np.unique(goodx) + b1, label="y={}x+{}".format(round(m1, 2), round(b1, 2)))
plt.plot(np.unique(badx), m2*np.unique(badx) + b2, label="y={}x+{}".format(round(m2, 2), round(b2, 2)))
#popt, pcov = curve_fit(func, np.asarray(x), np.asarray(y))
#plt.plot(x, func(np.asarray(x), *popt),label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
#plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.title("Surface Accessibility", fontsize=14, fontweight='bold')
plt.show()
fig.savefig("surface_RMSD.png") '''
################################################################
#GRAPH: Ca-RMSD vs charge
#DO FITTING!!!!!!!!!!!!!!!!
#################################################################

columns=["charge", "ID", "Surface", "Hydropathy", "local_AA","local_CA","global_AA","global_CA"]
data = pandas.read_csv(sys.argv[1], names=columns)
x=data.charge.tolist()
y=data.global_AA.tolist()


del x[0]
del y[0]
total=len(x)
x = [float(k) for k in x]
y = [float(k) for k in y]
gooddic={}
baddic={}

for i in range(len(x)):
	if y[i] > 2:
		if x[i] not in baddic:
			baddic.update({x[i]:[y[i]]})
		else:
			baddic[x[i]].append(y[i])
	else:
		if x[i] not in gooddic:
			gooddic.update({x[i]:[y[i]]})
		else:
			gooddic[x[i]].append(y[i])

dic_numbers={}
for keys in gooddic.keys():
	dic_numbers[keys]=len(gooddic[keys])
for keys in baddic.keys():
	if keys not in dic_numbers.keys():
			dic_numbers.update({keys:len(baddic[keys])})
	else:
		dic_numbers.update({keys: (dic_numbers[keys]+len(baddic[keys])) })

for keys in gooddic.keys():
	gooddic[keys]=len(gooddic[keys])/dic_numbers[keys]
for keys in baddic.keys():
	baddic[keys]=len(baddic[keys])/dic_numbers[keys]

listgood= (gooddic.items()) 
x_good, y_good= zip(*listgood) 
x_good=list(x_good)
y_good=list(y_good)
print(sum(y_good))

listbad= (baddic.items()) 
x_bad, y_bad= zip(*listbad) 
x_bad=list(x_bad)
y_bad=list(y_bad)
m1, b1 = np.polyfit(x_good, y_good, 1)
m2, b2 = np.polyfit(x_bad, y_bad, 1)
print(sum(y_bad))

#print(r1)
fig = plt.figure()
plt.xlabel('charge c')
plt.ylabel('%  of structures (with charge c)')
plt.plot(x_good, y_good,'o',markersize=5)
plt.plot(x_bad, y_bad,'ro',markersize=5)
plt.text(4,0.2, "RMSD <2", fontsize=10, fontweight='bold', color="blue")
plt.text(4,0.17, "RMSD >2", fontsize=10, fontweight='bold', color="red")
#plt.plot(np.unique(x_good), m1*np.unique(x_good) + b1, label="y={}x+{}".format(round(m1, 2), round(b1, 2)))
#plt.plot(np.unique(x_bad), m2*np.unique(x_bad) + b2, label="y={}x+{}".format(round(m2, 2), round(b2, 2)))
#popt, pcov = curve_fit(func, np.asarray(x), np.asarray(y))
#plt.plot(x, func(np.asarray(x), *popt),label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
#plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.title("Model Quality at Differing Charges", fontsize=14, fontweight='bold')
plt.show()
fig.savefig("charge_RMSD.png") 
'''

###############################################################
#GRAPH: Ca-RMSD vs charge
#DO FITTING!!!!!!!!!!!!!!!!
#################################################################

columns=["charge", "ID", "Surface", "Hydropathy", "local_AA","local_CA","global_AA","global_CA"]
data = pandas.read_csv(sys.argv[1], names=columns)
x=data.global_AA.tolist()

del x[0]

total=len(x)
x = [float(k) for k in x]
print(min(x), max(x))#0.24-25.1
steps=np.arange(0.2,25,0.3)
steps2=np.arange(0.5,25,0.3)
dic={}
for i in x:
	for a, b in zip(steps, steps2):
		if i>a and i<b:
			value=(a+b/2)
			if value in dic.keys():
				dic[value].append(i)
			else:
				dic.update({value: [i]})
		else:
			pass

for key, value in dic.items():
	dic[key]=len(dic[key])/total


list1 = (dic.items()) 
x, y= zip(*list1) 
x=list(x)
y=list(y)

fig = plt.figure()
plt.xlabel('RMSD global AA')
plt.ylabel('%  of structures')
plt.plot(x, y,'o',markersize=5)

#plt.text(4,0.15, "RMSD <2", fontsize=10, fontweight='bold', color="blue")
#plt.text(4,0.17, "RMSD >2", fontsize=10, fontweight='bold', color="red")
#plt.plot(np.unique(x_good), m1*np.unique(x_good) + b1, label="y={}x+{}".format(round(m1, 2), round(b1, 2)))
#plt.plot(np.unique(x_bad), m2*np.unique(x_bad) + b2, label="y={}x+{}".format(round(m2, 2), round(b2, 2)))
#popt, pcov = curve_fit(func, np.asarray(x), np.asarray(y))
#plt.plot(x, func(np.asarray(x), *popt),label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
#plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.title("RMSD Distribution", fontsize=14, fontweight='bold')
plt.show()
fig.savefig("distribution_RMSD.png") 
'''