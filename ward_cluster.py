import matplotlib.pyplot as plt
import pandas as pd
#%matplotlib inline
import numpy as np
from sklearn.cluster import AgglomerativeClustering
import sys
import scipy.cluster.hierarchy as shc
from scipy.cluster.hierarchy import ClusterWarning
from warnings import simplefilter
import math
import os
from scipy.spatial.distance import pdist, squareform
import itertools
from sklearn.cluster import DBSCAN
import collections
simplefilter("ignore", ClusterWarning)

def delete_nans(data):
	length=data.shape
	df2=data.replace(to_replace="NaN", value=np.nan, inplace=False)
	length_var=df2.shape[1]
	threshold=math.ceil(length_var*0.6)
	df2=df2.dropna(axis=1, how='all', thresh=threshold)
	df2=df2.dropna(axis=0, how='all', thresh=threshold)
	#print(df2)
	return df2
def get_dendograms(matrix_directory):
	list1 = os.listdir(matrix_directory) # dir is your directory path
	#print(len(list1))
	number_files = len(list1)
	if number_files % 2 != 0:
		number_files+=1
	plot_length=math.sqrt(number_files)
	plot_length=math.ceil(plot_length)
	plot_length=5
	counter=0
	shc.set_link_color_palette(['m', 'c', 'y', 'k'])
	fig, axes = plt.subplots(plot_length, plot_length, figsize=(20, 10))
	fig.tight_layout(pad=2.0)
	fig.subplots_adjust(top=0.88)
	#print("plot_length",plot_length)
	for a in range(plot_length):
		if counter<len(list1):
			for b in range(plot_length):
				drawn=0
				while drawn==0:
					if counter<len(list1):
						data = pd.read_csv(os.path.join(matrix_directory,list1[counter]), header=0, index_col=None)
						data=delete_nans(data)
						if data.shape[1]<=2:
							#print(a,b)
							#print(data)
							print(list1[counter])
						else:
							#print(a,b)
							data=pdist(data)
							dend = shc.dendrogram(shc.linkage(data, method='ward'),ax=axes[a,b])
							name=str(list1[counter]).replace("matrix_","")
							name=name.replace(".csv","")
							axes[a, b].set_title('Loop length n={}'.format(name))
							#axes[a,b].set(xlabel='x-label', ylabel='y-label')
							drawn=1
						counter+=1
						#print(counter)
					else:
						break
		else:break

	shc.set_link_color_palette(None)  # reset to default after use
	plt.suptitle("Structure Clustering of Same-Length Loops", size=16)
	plt.savefig("ward_dendogram.png")
	plt.show()

#get_dendograms(sys.argv[1])


def cluster_features(matrix_directory,combined_file):
	for file in os.listdir(os.path.join(matrix_directory)):
		clusterdic={"22":2, "12":3, "7":2, "18":3, "3":2, "6":2, "8":2, "21":2, "17":2, "10":2, "4":2, "16":3, "13":2, "9":2, "24":2, "15":2, "19":2, "5":2, "11":3,"20":2, "25":2, "14":3}
		data2 = pd.read_csv(os.path.join(combined_file), header=0, index_col=None)

		data = pd.read_csv(os.path.join(matrix_directory,file), header=0, index_col=None)
		data=delete_nans(data)
		#data=pdist(data)
		name=str(file).replace("matrix_","")
		name=name.replace(".csv","")
		identity=[]
		total_charge=[]
		nr_charged=[]
		Happiness_mean=[]
		Nr_sad=[]
		similarity=[]
		Hydropathy=[]
		Hydropathy_diff=[]
		Access=[]
		Relacc=[]
		Scacc=[]
		Screlacc=[]
		Access_avg=[]
		Relacc_avg=[]
		Scacc_avg=[]
		Screlacc_avg=[]
		simlength=[]

		feature_list=[identity,total_charge,nr_charged,Happiness_mean,Nr_sad,similarity,Hydropathy,Hydropathy_diff,Access,Relacc,Scacc,Screlacc,Access_avg,Relacc_avg,Scacc_avg,Screlacc_avg,simlength]
		feature_string=["identity","total_charge","nr_charged","Happiness_mean","Nr_sad","similarity","Hydropathy","Hydropathy_diff","Access","Relacc","Scacc","Screlacc","Access_avg","Relacc_avg","Scacc_avg","Screlacc_avg","simlength"]
		dic={"identity":identity,"total_charge":total_charge,"nr_charged":nr_charged,"Happiness_mean":Happiness_mean,"Nr_sad":Nr_sad,"similarity":similarity,"Hydropathy":Hydropathy,"Hydropathy_diff":Hydropathy_diff,"Access":Access,"Relacc":Relacc,"Scacc":Scacc,"Screlacc":Screlacc,"Access_avg":Access_avg,"Relacc_avg":Relacc_avg,"Scacc_avg":Scacc_avg,"Screlacc_avg":Screlacc_avg,"simlength":simlength}
		#dic=collections.OrderedDict()
		#print(feature_string)
		for i in data.columns:
			row=data2.loc[data2['ID'] == str(i)]
			identity+=row.identity.tolist()
			total_charge+=row.total_charge.tolist()
			nr_charged+=row.nr_charged.tolist()
			Happiness_mean+=row.Happiness_mean.tolist()
			Nr_sad+=row.Nr_sad.tolist()
			similarity+=row.similarity.tolist()
			Hydropathy+=row.Hydropathy.tolist()
			Hydropathy_diff+=row.Hydropathy_diff.tolist()
			Access+=row.Access.tolist()
			Relacc+=row.Relacc.tolist()
			Scacc+=row.Scacc.tolist()
			Screlacc+=row.Screlacc.tolist()
			Access_avg+=row.Access_avg.tolist()
			Relacc_avg+=row.Relacc_avg.tolist()
			Scacc_avg+=row.Scacc_avg.tolist()
			Screlacc_avg+=row.Screlacc_avg.tolist()
			simlength+=row.simlength.tolist()


		cluster = AgglomerativeClustering(n_clusters=clusterdic[name], affinity='euclidean', linkage='ward')
		cluster.fit_predict(data)
		print(data)
		print(cluster.labels_)
		#fig, axes = plt.subplots(len(feature_string), len(feature_string), figsize=(100, 100))
		#fig.tight_layout()
		#fig.subplots_adjust(top=0.88)
		combos=list(itertools.combinations(feature_string, 2))
		length_plt=(len(combos))
		length_plt=math.sqrt(length_plt)
		length_plt=math.ceil(length_plt)
		counter1=0
		counter2=0
		fig, axes = plt.subplots(length_plt, length_plt)
		fig.tight_layout()
		a=0
		b=0
		for y in range(0,length_plt):
			a=int(y)
			a=feature_string[a]
			for i in range(0,length_plt):
				b=int(i)
				b=feature_string[b]
				print(i,y)
				axes[i, y].scatter(dic[a],dic[b], c=cluster.labels_, cmap='rainbow')
				if i==(length_plt-1):
					axes[i,y].set(xlabel=a)
				if y==0:
					axes[i,y].set(ylabel=b)
				title=(str(feature_string[counter1])+str(feature_string[counter2]))
		plt.show()
		'''fig=plt.figure()
			print(a,b)
			print(dic[a],dic[b])
			plt.scatter(dic[a],dic[b], c=cluster.labels_, cmap='rainbow')
			plt.xlabel(a)
			plt.ylabel(b)
			plt.show()
		for a in range(len(feature_string)):
			if counter1!=0:
				counter2+=1
			counter1=0
			if counter1==maximum or counter2==maximum:
				break
			for b in range(len(feature_string)):
				if counter1==maximum or counter2==maximum:
					break
				axes[a, b].scatter(feature_list[counter1],feature_list[counter2], c=cluster.labels_, cmap='rainbow')
				axes[a,b].set(xlabel=feature_string[counter1], ylabel=feature_string[counter2])
				title=(str(feature_string[counter1])+str(feature_string[counter2]))
				print(feature_string[counter1],a,b)
				print(feature_string[counter2],a,b)

				#axes[a, b].set_title(str(feature_string[counter1])+str(feature_string[counter2]))
				#axes[a,b].axis('off')

				counter1+=1

		plt.show()'''

cluster_features(sys.argv[1],sys.argv[2])

def dbscan(matrix_directory):
	for file in os.listdir(os.path.join(matrix_directory)):
		clusterdic={"22":2, "12":3, "7":2, "18":3, "3":2, "6":2, "8":2, "21":2, "17":2, "10":2, "4":2, "16":3, "13":2, "9":2, "24":2, "15":2, "19":2, "5":2, "11":3,"20":2, "25":2, "14":3}
		data = pd.read_csv(os.path.join(matrix_directory,file), header=0, index_col=None)
		data=delete_nans(data)
		clustering = DBSCAN(eps=2, min_samples=2).fit(data)
		clustering.labels_
		print(clustering.labels_)

#dbscan(sys.argv[1])






#python3 ward_cluster.py ~/sync_project/Feature/length_matrix/ ~/sync_project/Feature/new.csv