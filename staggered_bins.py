#!/usr/bin/env python
import sys
import os
import pandas as pd
import math
from pathlib import Path
import numpy as np

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import style

# Evaluations
from sklearn.metrics import classification_report,confusion_matrix
# Random Forest
from sklearn.ensemble import RandomForestClassifier
# TensorFlow
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Flatten, Conv2D, MaxPooling2D
from tensorflow.keras.losses import sparse_categorical_crossentropy
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import KFold
from sklearn.metrics import matthews_corrcoef
from tensorflow.keras.models import Sequential, save_model
import joblib
from sklearn.preprocessing import MinMaxScaler
import collections
import shutil
from shutil import copyfile
'''original setting:
	staggering_value=0.1
	max_angstrom=5
	first_layer_size=0.5
	second_layer_size=1
	'''

def RMSD_binary(threshold,RMSD_df):
	columns=["local_AA_bin"+str(threshold),"local_CA_bin"+str(threshold),"global_AA_bin"+str(threshold), "global_CA_bin"+str(threshold), "ID"]
	data=[]
	local_AA= RMSD_df.local_AA.tolist()
	local_CA= RMSD_df.local_CA.tolist()
	global_AA= RMSD_df.global_AA.tolist()
	global_CA=RMSD_df.global_CA.tolist()
	name=RMSD_df.ID.tolist()
	for (a,b,c,d,n) in zip(local_AA, local_CA, global_AA, global_CA, name):
		a=float(a)
		b=float(b)
		c=float(c)
		d=float(d)
		if a>=threshold:
			a=1
		else:
			a=0
		if b>=threshold:
			b=1
		else:
			b=0
		if c>=threshold:
			c=1
		else:
			c=0
		if d>=threshold:
			d=1
		else:
			d=0
		write=[a,b,c,d,n]
		data.append(write)
	dfnew = pd.DataFrame(data, columns=columns)
	#print(dfnew)
	return (dfnew)

def RMSD_nom(RMSD_df,second_layer_list_up,second_layer_list_down):
	columns=["local_AA_nom","local_CA_nom","global_AA_nom", "global_CA_nom", "ID"]
	data=[]
	local_AA= RMSD_df.local_AA.tolist()
	local_CA= RMSD_df.local_CA.tolist()
	global_AA= RMSD_df.global_AA.tolist()
	global_CA=RMSD_df.global_CA.tolist()
	name=RMSD_df.ID.tolist()
	for (a,b,c,d,n) in zip(local_AA, local_CA, global_AA, global_CA, name):
		for (threshold_up,threshold_low) in zip(second_layer_list_up, second_layer_list_down):
			#print(threshold_low,threshold_up)
			#print(type(a))
			if threshold_up==max(second_layer_list_up):
				threshold_up=50
			if type(a)==str:
				pass
			elif a>=threshold_low and a<threshold_up:
				a=">"+str(threshold_low)+"<"+str(threshold_up)
			if type(b)==str:
				pass
			elif b>=threshold_low and b<threshold_up:
				b=">"+str(threshold_low)+"<"+str(threshold_up)
			if type(c)==str:
				pass
			elif c>=threshold_low and c<threshold_up:
				c=">"+str(threshold_low)+"<"+str(threshold_up)
			if type(d)==str:
				pass
			elif d>=threshold_low and d<threshold_up:
				d=">"+str(threshold_low)+"<"+str(threshold_up)
		write=[a,b,c,d,n]
		#print(write)
		data.append(write)
	dfnew = pd.DataFrame(data, columns=columns)
	#print(dfnew)
	return (dfnew)

def make_bin_df(RMSD_file, feature_csv, bin_models,staggering_value,max_angstrom,first_layer_size,second_layer_size,accuracy_csv_path):
	#here the staggering is set
	#total range: 0-5 i.e above 5Angstrom there is no more differentiation
	#first_layer bins are 1 Angstrom in size, second layer bins are 2Angstroms in size
	#staggering value determines the shift from bin to bin in the first layer
	second_layer_list_down=np.arange(0,max_angstrom,second_layer_size) # makes 0,1,2,...4
	second_layer_list_up=np.arange(second_layer_size,max_angstrom,second_layer_size) # makes 1,2,...4
	second_layer_list_up=np.append(second_layer_list_up,max_angstrom)
	first_layer_thresh_list=[]
	RMSD_path=os.path.abspath(RMSD_file)
	feature_path=os.path.abspath(feature_csv)
	RMSD_df=pd.read_csv(RMSD_path,header=0)
	feature_df=pd.read_csv(feature_path,header=0)
	print("__________________________noms:up, down", list(second_layer_list_up),list(second_layer_list_down))

	df=RMSD_nom(RMSD_df,list(second_layer_list_up),list(second_layer_list_down))
	for (up,down) in zip(second_layer_list_up,second_layer_list_down):
		staggerbins_down=np.arange(down,up,staggering_value)
		bin_list=[]
		x=first_layer_size
		for i in staggerbins_down:
			if x<up:
				x=i+first_layer_size
				x=np.around(x,2)
				i=np.around(i,2)
				bin_list.append((i,x))
		first_layer_thresh_list.append(bin_list)
	first_layer_thresh_list_down=np.arange(0,(max_angstrom-first_layer_size), staggering_value )
	first_layer_thresh_list_down=[np.around(x,2) for x in first_layer_thresh_list_down]
	first_layer_thresh_list_up=np.arange(first_layer_size, max_angstrom, staggering_value)
	first_layer_thresh_list_up=[np.around(x,2) for x in first_layer_thresh_list_up]
	first_layer_thresh_list=zip(first_layer_thresh_list_down,first_layer_thresh_list_up)
	print(first_layer_thresh_list)
	accuracy_bins={}
	for down,up in first_layer_thresh_list:
		print("_____________________ bins",down,up)
		df2=RMSD_binary(up,RMSD_df)
		prediction_df=df2.merge(feature_df, on="ID")
		print(prediction_df)
		print(prediction_df.columns)
		prediction_df2=parse_prediction_df(prediction_df,up, bin_models)
		vals,bin_vals= Kfolds(prediction_df2,up, bin_models)
		mean_accuracy= model(vals, bin_vals, bin_models, up)
		accuracy_bins.update({up:[mean_accuracy]})
	accuracy_df = pd.DataFrame.from_dict(accuracy_bins)
	accuracy_df.to_csv(os.path.join(accuracy_csv_path), index=False)
	return df
	

def parse_prediction_df(prediction_df, up, bin_models):
	features_to_include=["local_CA_bin"+str(up),'LLE1', 'LLE2', 'LLE3', 'LLE4', 'LLE5', 'LLE6', 'LTSA1',
       'LTSA2', 'LTSA3', 'LTSA4', 'LTSA5', 'LTSA6', 'Mod_LLE1', 'Mod_LLE2',
       'Mod_LLE3', 'Mod_LLE4', 'Mod_LLE5', 'Mod_LLE6', 'Isomap1', 'Isomap2',
       'Isomap3', 'Isomap4', 'Isomap5', 'Isomap6', 'MDS1', 'MDS2', 'MDS3',
       'MDS4', 'MDS5', 'MDS6', 'SE1', 'SE2', 'SE3', 'SE4', 'SE5', 'SE6',
       'tSNE1', 'tSNE2', 'tSNE3', 'tSNE4', 'tSNE5', 'tSNE6', 'blosum62_val1',
       'blosum62_val2', 'blosum62_val3', 'blosum62_val4', 'blosum62_val5',
       'blosum62_val6', 'tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'Happiness_mean',
       'Nr_sad', 'identity', 'similarity','Hydropathy',
       'Hydropathy_diff', 'Access', 'Relacc', 'Scacc', 'Screlacc','simlength'] #deleted tip_res
	prediction_df.isnull().any()

	prediction_df2 = prediction_df.copy()
	prediction_df2 = pd.DataFrame(data=prediction_df2, columns=features_to_include)
	print("hey")
	print(prediction_df2.head())
	return prediction_df2
	

def Kfolds(prediction_df2,up, bin_models):
	bin_vals="local_CA_bin"+str(up)
	#input(prediction_df2.iloc[:,1:])
	vals=prediction_df2.iloc[:,1:].to_numpy()
	bin_vals= prediction_df2[[bin_vals]].to_numpy()
	#bin_vals=tuple(bin_vals)
	#vals=tuple(vals)
	#input(bin_vals)
	#input(vals)
	return vals,bin_vals

def model(X,y, bin_models,up):
	num_folds = 10
	acc_per_fold = []
	MCC_per_fold=[]
	fold_no = 1
	kfold = KFold(n_splits=num_folds, shuffle=True)
	fold_dic={}
	os.mkdir(os.path.join(bin_models,"folds"))
	for train, test in kfold.split(X, y):
		print('------------------------------------------------------------------------')
		print(f'Training for fold {fold_no} ...')
		model = RandomForestClassifier(n_estimators=200)
		model.fit(X[train], y[train])
		RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
			max_depth=None, max_features='auto', max_leaf_nodes=None,
			min_impurity_decrease=0.0, min_impurity_split=None,
			min_samples_leaf=1, min_samples_split=2,
			min_weight_fraction_leaf=0.0, n_estimators=200, n_jobs=None,
			oob_score=False, random_state=None, verbose=0,
			warm_start=False)

		rf_prediction = model.predict(X[test])# Evaluations
		print('Classification Report: \n')
		MCC_of_fold=matthews_corrcoef(y[test],rf_prediction, sample_weight=None )
		print("MCC",MCC_of_fold)

		print(classification_report(y[test],rf_prediction))
		print('\nConfusion Matrix: \n')
		print(confusion_matrix(y[test],rf_prediction))
		#input("complete")
		#########################
		scores = model.score(X[test], y[test])
		print(scores)
		#print(f'Score for fold {fold_no}: {model.metrics_names[0]} of {scores[0]}; {model.metrics_names[1]} of {scores[1]*100}%')
		acc_per_fold.append(scores * 100)
		MCC_per_fold.append(MCC_of_fold)
		joblib.dump(model, os.path.join(bin_models,"folds","bin"+str(up)+"fold"+str(fold_no)+".pkl")) 
		fold_dic.update({MCC_of_fold:fold_no})
		fold_no = fold_no + 1
	max_MCC=max(MCC_per_fold)
	max_fold=fold_dic[max_MCC]
	copyfile(os.path.join(bin_models,"folds","bin"+str(up)+"fold"+str(max_fold)+".pkl"), os.path.join(bin_models, "bin"+str(up)+".pkl"))
	shutil.rmtree(os.path.join(bin_models,"folds"), ignore_errors=True)
	#joblib.dump(model, os.path.join(bin_models, "bin"+str(up)+".pkl")) 

	# == Provide average scores ==
	print('------------------------------------------------------------------------')
	print('Score per fold')
	for i in range(0, len(acc_per_fold)):
	  print('------------------------------------------------------------------------')
	  print(f'> Fold {i+1} -  Accuracy: {acc_per_fold[i]}%')
	print('------------------------------------------------------------------------')
	print('Average scores for all folds:')
	print(f'> Accuracy: {np.mean(acc_per_fold)} (+- {np.std(acc_per_fold)})')
	print('------------------------------------------------------------------------')


	return np.mean(MCC_per_fold)

def MCC_weighting(accuracy_bins,second_layer_list_up, predicted_vals, first_layer_size,staggering_value,max_angstrom):
	accuracies=collections.OrderedDict()
	scaler = MinMaxScaler()                   ##############dont do min max scaling!!! use 0 as min and max as max
	for i in second_layer_list_up:
		i_low=i-1
		accuracy_of_key=[]
		for key,value in predicted_vals.items():
			key=float(key)
			#print(accuracy_bins)
			if key <=i and key>i_low:
				accuracy_of_key.append(accuracy_bins[str(key)][0])
		accuracies.update({i:accuracy_of_key})
	print(accuracies)
	for key,value in accuracies.items():
		value=np.array(value)
		value=value.reshape(-1, 1)
		print("values",value)
		scaler.fit(value)
		value=scaler.transform(value)
		accuracies[key]=value  
	accuracy_list=[]
	for key,value in accuracies.items():
		value=value.tolist()
		accuracy_list+=value
	first_layer_thresh_list_up=np.arange(first_layer_size, max_angstrom, staggering_value)
	first_layer_thresh_list_up=[np.around(x,2) for x in first_layer_thresh_list_up]
	new_list=zip(first_layer_thresh_list_up,accuracy_list)
	new_dic={}
	print(first_layer_thresh_list_up)
	print(accuracy_list)
	print(new_list)
	for (key,value) in new_list:
		print(key,value)
		new_dic.update({key:value})
	print(new_dic)
	return new_dic


def make_test_df(bin_models, test_file,second_layer_size, max_angstrom,first_layer_size,staggering_value): #accuracy_csv_path
	features_to_include=['LLE1', 'LLE2', 'LLE3', 'LLE4', 'LLE5', 'LLE6', 'LTSA1',
       'LTSA2', 'LTSA3', 'LTSA4', 'LTSA5', 'LTSA6', 'Mod_LLE1', 'Mod_LLE2',
       'Mod_LLE3', 'Mod_LLE4', 'Mod_LLE5', 'Mod_LLE6', 'Isomap1', 'Isomap2',
       'Isomap3', 'Isomap4', 'Isomap5', 'Isomap6', 'MDS1', 'MDS2', 'MDS3',
       'MDS4', 'MDS5', 'MDS6', 'SE1', 'SE2', 'SE3', 'SE4', 'SE5', 'SE6',
       'tSNE1', 'tSNE2', 'tSNE3', 'tSNE4', 'tSNE5', 'tSNE6', 'blosum62_val1',
       'blosum62_val2', 'blosum62_val3', 'blosum62_val4', 'blosum62_val5',
       'blosum62_val6', 'tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'Happiness_mean',
       'Nr_sad', 'identity', 'similarity','Hydropathy',
       'Hydropathy_diff', 'Access', 'Relacc', 'Scacc', 'Screlacc','simlength'] #deleted tip_res
	test_df=pd.read_csv(test_file)
	test_df.isnull().any()
	#accuracy_bins=pd.read_csv(accuracy_csv_path, header=0)
	#accuracy_bins=accuracy_bins.to_dict()

	test_df2 = test_df.copy()
	test_df2 = pd.DataFrame(data=test_df2, columns=features_to_include)
	#print(test_df2.head())
	predicted_vals={}
	for file in os.listdir(bin_models):
		filename=file.replace(".pkl","")
		filename=filename.replace("bin","")
		model= joblib.load(os.path.join(bin_models,file))
		vals=test_df2.iloc[:,:].to_numpy()
		predictions=model.predict(vals)
		predicted_vals.update({filename:predictions})
	
	second_layer_list_up=np.arange(second_layer_size,max_angstrom,second_layer_size) # makes 1,2,...4
	second_layer_list_up=np.append(second_layer_list_up,5)
	summed_layer_two={}
	#accuracy_dic=MCC_weighting(accuracy_bins,second_layer_list_up, predicted_vals, first_layer_size,staggering_value,max_angstrom)
	for i in second_layer_list_up:
		i_low=i-1
		multiply_list=[]
		accuracies=[]
		for key,value in predicted_vals.items():
			key=float(key)
			if key <=i and key>i_low:
				#multiply_by=accuracy_dic[key]
				multiply_list+=[value] #*multiply_by
				#print("HEEEEEEEEEEEEEY",type(value))
		number_of_vals=len(multiply_list)
		arr = np.array(multiply_list)
		sum_array=arr.sum(axis=0)/number_of_vals   
		#print("sum array",sum_array)
		#sum_array = sum_array/number_of_vals
		summed_layer_two.update({i:sum_array})
	test_df=pd.DataFrame.from_dict(summed_layer_two)
	print("test_df",test_df)
	test_df.columns = test_df.columns.astype(str)
	return test_df



def third_layer(a, third_model_path):
	a.isnull().any()
	a_nom=a[["local_CA_nom"]].to_numpy()
	a_vals=a.drop(labels="local_CA_nom",axis=1).to_numpy()
	up="THIRD_LAYER"
	model(a_vals, a_nom, third_model_path, up)


def run_with_unknown(unknown_file, third_model_path) :
	staggering_value=0.2
	max_angstrom=5
	first_layer_size=1
	second_layer_size=1
	features_to_include=['LLE1', 'LLE2', 'LLE3', 'LLE4', 'LLE5', 'LLE6', 'LTSA1',
       'LTSA2', 'LTSA3', 'LTSA4', 'LTSA5', 'LTSA6', 'Mod_LLE1', 'Mod_LLE2',
       'Mod_LLE3', 'Mod_LLE4', 'Mod_LLE5', 'Mod_LLE6', 'Isomap1', 'Isomap2',
       'Isomap3', 'Isomap4', 'Isomap5', 'Isomap6', 'MDS1', 'MDS2', 'MDS3',
       'MDS4', 'MDS5', 'MDS6', 'SE1', 'SE2', 'SE3', 'SE4', 'SE5', 'SE6',
       'tSNE1', 'tSNE2', 'tSNE3', 'tSNE4', 'tSNE5', 'tSNE6', 'blosum62_val1',
       'blosum62_val2', 'blosum62_val3', 'blosum62_val4', 'blosum62_val5',
       'blosum62_val6', 'ID', 'tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'Happiness_mean',
       'Nr_sad', 'identity', 'similarity', 'Hydropathy',
       'Hydropathy_diff', 'Access', 'Relacc', 'Scacc', 'Screlacc',
       'local_CA','simlength'] #deleted tip_res
	second_layer_list_down=np.arange(0,max_angstrom,second_layer_size) # makes 0,1,2,...4
	second_layer_list_up=np.arange(second_layer_size,max_angstrom,second_layer_size) # makes 1,2,...4
	second_layer_list_up=np.append(second_layer_list_up,max_angstrom)
	bin_models=Path("/home/lilian/sync_project/staggered_models")
	#first we take our csv with all our features and remove all unwanted features (and calculate local_CA_nom class if you want for testing)

	#we then feed into first layer and calculate the prediction values for each categories from second layer back into feature csv
	prediction_df=make_test_df(bin_models, unknown_file, second_layer_size,max_angstrom,first_layer_size,staggering_value)
	unknown_file=pd.read_csv(unknown_file,header=0)
	df_nom=RMSD_nom(unknown_file,second_layer_list_up,second_layer_list_down)
	df_nom=df_nom[["local_CA_nom","ID"]]
	#using third layer model we make final prediction
	unknown_file=unknown_file.merge(df_nom,on="ID")
	a=prediction_df.join(unknown_file,on=None)
	print(a.columns)
	second_layer_list_up=np.arange(second_layer_size,max_angstrom,second_layer_size)
	second_layer_list_up=[str(x) for x in second_layer_list_up]
	second_layer_list_up.append(str(max_angstrom))
	print(second_layer_list_up)
	list_features=second_layer_list_up+features_to_include+["local_CA_nom"]
	a=a[list_features]
	#a=a.merge(nom_df,on="ID")
	a=a.drop(labels=["ID"], axis=1)
	a.to_csv(os.path.join("/home/lilian/sync_project/Feature","train_staggered_pre.csv"), index=False)
	print(a.head())
	model= joblib.load(os.path.join(third_model_path))
	nom=a[["local_CA_nom"]].to_numpy()
	vals=a.drop(labels="local_CA_nom",axis=1).to_numpy()
	predictions=model.predict(vals)

	print('Classification Report: \n')
	MCC_of_fold=matthews_corrcoef(nom,predictions, sample_weight=None )
	print("MCC",MCC_of_fold)

	print(classification_report(nom,predictions))
	print('\nConfusion Matrix: \n')
	print(confusion_matrix(nom,predictions))
	#input("complete")
	#########################
	scores = model.score(vals, nom)
	print(scores)

	plt.figure(figsize=(12,5))
	sns.heatmap(a.corr(),annot=True)
	plt.show()
	plt.figure(figsize=(12,5))
	sns.countplot(x='local_CA_nom', data=a)
	plt.show()
	plt.figure(figsize=(12,5))
	sns.pairplot(data=a, hue='local_CA_nom')
	plt.show()




def run(RMSD_file, feature_csv, bin_models, test_file,feature_directory,accuracy_csv_path,third_model_path):
	staggering_value=0.2
	max_angstrom=5
	first_layer_size=1
	second_layer_size=1
	features_to_include=['LLE1', 'LLE2', 'LLE3', 'LLE4', 'LLE5', 'LLE6', 'LTSA1',
       'LTSA2', 'LTSA3', 'LTSA4', 'LTSA5', 'LTSA6', 'Mod_LLE1', 'Mod_LLE2',
       'Mod_LLE3', 'Mod_LLE4', 'Mod_LLE5', 'Mod_LLE6', 'Isomap1', 'Isomap2',
       'Isomap3', 'Isomap4', 'Isomap5', 'Isomap6', 'MDS1', 'MDS2', 'MDS3',
       'MDS4', 'MDS5', 'MDS6', 'SE1', 'SE2', 'SE3', 'SE4', 'SE5', 'SE6',
       'tSNE1', 'tSNE2', 'tSNE3', 'tSNE4', 'tSNE5', 'tSNE6', 'blosum62_val1',
       'blosum62_val2', 'blosum62_val3', 'blosum62_val4', 'blosum62_val5',
       'blosum62_val6', 'ID', 'tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'Happiness_mean',
       'Nr_sad', 'identity', 'similarity', 'Hydropathy',
       'Hydropathy_diff', 'Access', 'Relacc', 'Scacc', 'Screlacc',
       'local_CA','simlength'] #deleted tip_res

	nom_df=make_bin_df(RMSD_file, feature_csv, bin_models,staggering_value,max_angstrom,first_layer_size,second_layer_size, accuracy_csv_path)
	test_df=make_test_df(bin_models, test_file, second_layer_size,max_angstrom,first_layer_size,staggering_value)
	feature_csv=pd.read_csv(os.path.abspath(feature_csv))
	print(feature_csv)
	a=test_df.join(feature_csv,on=None)
	print(a.columns)
	second_layer_list_up=np.arange(second_layer_size,max_angstrom,second_layer_size)
	second_layer_list_up=[str(x) for x in second_layer_list_up]
	second_layer_list_up.append(str(max_angstrom))
	print(second_layer_list_up)
	list_features=second_layer_list_up+features_to_include
	a=a[list_features]
	a=a.merge(nom_df,on="ID")
	a=a.drop(labels=["ID","local_AA_nom","global_AA_nom","global_CA_nom"], axis=1)
	print(a.head())
	third_layer(a, third_model_path)
	a.to_csv(os.path.join(feature_directory,"MERGED_STAGGERED_tmp"), index=False)








run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[5])
#run_with_unknown(sys.argv[1],sys.argv[2])
#run_with_unknown(sys.argv[3],sys.argv[4])



#with setting:

'''staggering_value=0.5
	max_angstrom=5
	first_layer_size=1
	second_layer_size=2'''
#####=========0.6266835110483702

'''staggering_value=0.2
	max_angstrom=5
	first_layer_size=1
	second_layer_size=2'''
#####=========MCC=0.6471775637524335 		!!!!

'''staggering_value=0.2
	max_angstrom=5
	first_layer_size=1.5
	second_layer_size=3'''
####=========0.6274200759280766
######################################################################################

'''staggering_value=0.2
	max_angstrom=5
	first_layer_size=1.5
	second_layer_size=2'''
####MCC=0.6572073279745841


