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
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
#################################################################
#GRAPH: Sequence Similarity/Identity vs RMSD
#DO FITTING!!!!!!!!!!!!!!!!
#################################################################
data = pd.read_csv(sys.argv[1])
x=data.identity.tolist()
y=data.simlength.tolist()
a=data.local_CA.tolist()
length=data.length.tolist()
length=[int(x) for x in length]

scaler = MinMaxScaler() 
column_names_to_normalize = ['simlength']
x = data[column_names_to_normalize].values
x_scaled = scaler.fit_transform(x)
data_temp = pd.DataFrame(x_scaled, columns=column_names_to_normalize, index = data.index)
data[column_names_to_normalize] = data_temp

data = data[data["simlength"].between(0.3, 100)] #get rid of the two outliers
data['local_CA_log'] = np.log10(data['local_CA'])


df1 = pd.DataFrame(list(zip(x, a)), columns =['sequence identity', 'RMSD Local CA (Å)']) 
df2 = pd.DataFrame(list(zip(y, a)), columns =['sequence similarity', 'RMSD Local CA (Å)']) 
sns.jointplot(data=data, x="identity", y="local_CA", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}},scatter_kws={'alpha':0.0})
ax=sns.scatterplot(data=data, x="identity", y="local_CA",hue="length", palette='viridis',alpha=0.5, s=50)
sm = plt.cm.ScalarMappable(cmap="viridis")
sm.set_array([])
plt.xlim([-0.05,1.01])
plt.ylim([-0.1,8.5])
ax.get_legend().remove()
ticks=[0,2,4,6,8,10, 12,14,16,18,20,22,24,26,28]
nticks=[((x - min(length)) / (max(length) - min(length))) for x in ticks]
ticks=[str(x) for x in ticks]
#cbar=ax.figure.colorbar(sm, ticks=nticks, orientation='horizontal')
#cbar.ax.set_xticklabels(ticks) 
plt.savefig("identity_color.png")

sns.jointplot(data=data, x="simlength", y="local_CA_log", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}},scatter_kws={'alpha':0.0})
ax=sns.scatterplot(data=data, x="simlength", y="local_CA_log",hue="length", palette='viridis',alpha=0.5, s=50)
sm = plt.cm.ScalarMappable(cmap="viridis")
sm.set_array([])
#plt.xlim([-0.1,0.3])
#plt.ylim([0,8.5])
ax.get_legend().remove()
ticks=[0,2,4,6,8,10, 12,14,16,18,20,22,24,26,28]
nticks=[((x - min(length)) / (max(length) - min(length))) for x in ticks]
ticks=[str(x) for x in ticks]
#cbar=ax.figure.colorbar(sm, ticks=nticks)
#cbar.ax.set_yticklabels(ticks) 
plt.savefig("similarity_length_color.png")

sns.jointplot(data=data, x="identity", y="local_CA", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
plt.xlim([-0.05,1.01])
plt.ylim([-0.1,8.5])
plt.savefig("identity.png") 

sns.jointplot(data=data, x="simlength", y="local_CA_log", kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}})
#plt.xlim([-0.6,1.01])
#plt.ylim([0,8.5])
plt.savefig("similarity_length.png") 

plt.show()