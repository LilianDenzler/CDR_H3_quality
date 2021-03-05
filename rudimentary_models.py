#!/usr/bin/env python
import sys
import os
import pandas as pd
import math
from pathlib import Path
import numpy as np

import pandas as pd
import qualiloop
from qualiloop import myfunctions
sys.path.append('./libs')
sys.path.append('./CDRH3lib')


##############################################################################
#Make Dataframe of RMSD values (numerical and nominal) and length, identity and similarity as features
def rudimentary_features(file, log_name):
	#get loop length
	loop_seq_df,length_seq_df=myfunctions.get_loop_seq (file)
	length_df=length_seq_df[["length"]]
	length_df=length_df.copy()
	#get similarity and identity
	identity_similarity_df=handle_log_file(log_name)
	id_sim_df=identity_similarity_df[["identity", "similarity"]]
	id_sim_df=id_sim_df.copy()
	rud_features_df=pd.concat([length_df, id_sim_df])
	return rud_features_df

def make_df_num(file, actual_file, threshold_list=[1,2,3,4,5,6,7,8,9,10], threshlow=[0,2,4,6,8,10,200], threshup=[2,4,6,8,10,12,200]):
	RMSD_num_df=myfunctions.get_RMSD_num(file, actual_file)
	df_all_bins=myfunctions.get_all_bins(file,actual_file, threshold_list)
	df_nom=myfunctions.RMSD_nom(threshlow,threshup, file, actual_file)
	#combine to one dataframe
	all_RMSD=pd.concat([RMSD_num_df, df_all_bins, df_nom])
	return all_RMSD

def get_full_dataset_df(model_dir, log_dir,actual_dir, threshold_list=[1,2,3,4,5,6,7,8,9,10], threshlow=[0,2,4,6,8,10,200], threshup=[2,4,6,8,10,12,200]):
	full_dataset_df=pd.Dataframe()
	for file in os.listdir(model_dir):
		filename=os.path.splitext(os.path.basename(file))[0]
		actual_file=os.path.join(actual_dir, filename+".pdb")
		log_name=os.path.join(log_dir, filename+".log")
		rud_features_df=rudimentary_features(file, log_name)
		all_RMSD=make_df_num(file, actual_file, threshold_list, threshlow, threshup)
		final_df=pd.concat([rud_features_df, all_RMSD])
		if full_dataset_df.empty != True:
			full_dataset_df=pd.concat([full_dataset_df, final_df])
		else:
			full_dataset_df=final_df
	return full_dataset_df

if __name__ == '__main__':
	model_dir=sys.argv[1]
	log_dir=sys.argv[2]
	actual_dir=sys.argv[3]
	threshold_list=[1,2,3,4,5,6,7,8,9,10]
	threshlow=[0,2,4,6,8,10,200]
	threshup=[2,4,6,8,10,12,200]
	full_dataset_df=get_full_dataset_df(model_dir, log_dir,actual_dir, threshold_list, threshlow, threshup)


		

		

	

