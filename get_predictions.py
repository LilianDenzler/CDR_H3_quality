import sys
import os
import pandas as pd
#import save_RMS_lib
#import get_abYmod2_lib
#import seq_parser_lib
from pathlib import Path
import numpy as np
import epitopepredict as ep
from skimage.feature import hog
from skimage import data, exposure
from sklearn.decomposition import TruncatedSVD
import collections
from sklearn import manifold, datasets
from functools import partial

def get_results(results_directory,featurefile):
	featurefile=pd.read_csv(featurefile,header=0)
	#feature_file=pd.read_csv(feature_file,header=0)
	#featurefile=featurefile[["local_CA_nom"]]
	print(featurefile)
	df_dict={}
	for file in os.listdir(results_directory):
		print(file)
		filename=file.replace("bin","")
		filename=file.replace("_results","")
		df=pd.read_csv(os.path.join(results_directory,file), header=0)
		print(df)
		df[['weka_nr_pred',"predicted"]] = df.predicted.str.split(":",expand=True)
		#df[['weka_nr_act',"actual"]] = df.actual.str.split(":",expand=True)
		df=df[["predicted", "prediction"]]
		df=df.rename(columns={"predicted": "predicted"+filename,"prediction":"prediction"+filename})
		print(df)
		df_dict.update({filename:df})
	keys=list(df_dict.keys())
	a=keys[0]
	a=df_dict[a]
	for i in keys:
		if keys.index(i)==0:
			continue
		b=df_dict[i]
		merged=pd.concat([a,b],axis=1)
		a=merged
	a=pd.concat([a,featurefile],axis=1)
	print(a)
	a.to_csv(os.path.join(results_directory,"feature+predictions"+".csv"), index=False)
	#merged=a.merge(feature_file, on="local_CA_nom")
	#print(merged)
	#merged.to_csv(os.path.join(results_directory,"full_feature"+".csv"), index=False)

get_results(sys.argv[1],sys.argv[2])
