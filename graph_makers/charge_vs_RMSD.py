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

data = pandas.read_csv(sys.argv[1])
x=data.total_charge.tolist()
a=data.nr_charged.tolist()
y=data.local_AA.tolist()

#Red: y/ (percentage of bad structures overall)
#blue: y/ (percentage of good structures overall)


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
print(dic_numbers)
for keys in baddic.keys():
	if keys not in dic_numbers.keys():
			dic_numbers.update({keys:len(baddic[keys])})
	else:
		dic_numbers.update({keys: (dic_numbers[keys]+len(baddic[keys])) })

totalgood=0
for keys in gooddic.keys():
	totalgood+=len(gooddic[keys])
print(totalgood)

totalbad=0
for keys in baddic.keys():
	totalbad+=len(baddic[keys])
print(totalbad)

total=0
for keys in dic_numbers.keys():
	total+=dic_numbers[keys]
print(totalbad)


for keys in gooddic.keys():
	print("percentage of good at charge", len(gooddic[keys])/dic_numbers[keys])
	print("percentage of total good", (totalgood/total))
	gooddic[keys]=(len(gooddic[keys])/dic_numbers[keys])/(totalgood/total)
	print(gooddic[keys])
for keys in baddic.keys():
	print("percentage of bad at charge", len(baddic[keys])/dic_numbers[keys])
	print("percentage of total bad", (totalbad/total))
	baddic[keys]=(len(baddic[keys])/dic_numbers[keys])/(totalbad/total)
	print(baddic[keys])

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
plt.ylabel('%  of structures (with charge c) ')
plt.plot(x_good, y_good,'o',markersize=5)
plt.plot(x_bad, y_bad,'ro',markersize=5)
plt.text(3,0.3, "RMSD <2", fontsize=10, fontweight='bold', color="blue")
plt.text(3,0.2, "RMSD >2", fontsize=10, fontweight='bold', color="red")
#plt.plot(np.unique(x_good), m1*np.unique(x_good) + b1, label="y={}x+{}".format(round(m1, 2), round(b1, 2)))
#plt.plot(np.unique(x_bad), m2*np.unique(x_bad) + b2, label="y={}x+{}".format(round(m2, 2), round(b2, 2)))
#popt, pcov = curve_fit(func, np.asarray(x), np.asarray(y))
#plt.plot(x, func(np.asarray(x), *popt),label='fit: a=%5.3f, b=%5.3f' % tuple(popt))
#plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.title("Model Quality at Differing Charges", fontsize=14, fontweight='bold')
plt.show()
fig.savefig("charge_RMSD.png") 