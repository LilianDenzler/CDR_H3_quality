#!/usr/bin/env python
import sys
import os
import numpy as np
import csv
import pandas as pd
import itertools



def sort_by_length(combined_csv, feature_directory):
	df=pd.read_csv(combined_csv,header=0)
	columns=df.columns
	columns=columns.tolist()
	length=df.length.tolist()
	length=[int(x) for x in length]
	max_length=max(length)
	min_length=min(length)
	dic={}
	for i in range(min_length,(max_length+1)):
		for index, row in df.iterrows():
			if row['length']==i:
				if str(i) in dic:
					dic[str(i)].append(row.tolist())
				else:
					dic[str(i)]=[row.tolist()]


	if os.path.exists(os.path.join(feature_directory,"length_features")):
		pass
	else:
		os.mkdir(os.path.join(feature_directory,"length_features"))

	for key in dic:
		name="length_combined_"+str(key)
		list_of_dic=dic[key]
		print(list_of_dic)
		df = pd.DataFrame(list_of_dic,columns=columns, index=None)
		df.to_csv(os.path.join(feature_directory,"length_features",name+'.csv'), index=None)

def get_RMSDs(feature_directory,script_name,actual_directory):

	for file in os.listdir(os.path.join(feature_directory,"length_features")):
		df=pd.read_csv(os.path.join(feature_directory,"length_features",file),header=0)
		names_id=df.ID.tolist()
		filename=file.replace("length_combined_", "")
		filename=filename.replace(".csv", "")
		for (a,b) in itertools.combinations(names_id, 2):

			command=os.popen("./profit -f {} -h {} {}".format(script_name, os.path.join(actual_directory,a+".pdb"), os.path.join(actual_directory,b+".pdb"))).readlines()
			if not os.path.exists(os.path.join(feature_directory,"length_RMSD_ByRes")):
				os.makedirs(os.path.join(feature_directory,"length_RMSD_ByRes"))
			By_res=open(os.path.join(feature_directory,"length_RMSD_ByRes",filename+ "_"+a+"_"+b+".profit_by_res"),'w+')
			counter=0
			list_write=[]
			for line in command:
				print(line)
				columns=line.split("    ")
				if len(columns)< 2 and 'RMS:' in line:
					line=line.replace(" ","")
					line=line.replace("\n","")
					list_write.append(line[4:])
				if len(columns)>=3 and 'RMS:' in line:
					By_res.write(line)
					counter=1
				elif counter==1:
					By_res.write(" \n")
					counter=0
			if len(list_write)>7:
				del list_write[2]
				del list_write[3]
				del list_write[3]
				del list_write[4]

			list_write.append(a)
			list_write.append(b)
			print(list_write)
			columns=["local_AA","local_CA","global_AA","global_CA","one","two"]
			if os.path.exists(os.path.join(feature_directory,"RMS_" +file+".csv")):
				with open(os.path.join(feature_directory,"RMS_" + file+".csv"), "a") as f:
					writer = csv.writer(f)
					writer.writerow(list_write)
			else:
				with open(os.path.join(feature_directory,"RMS_" + file+".csv"),"w") as f:
					writer = csv.writer(f)
					writer.writerow(columns)
					writer.writerow(list_write)

def get_matrix(feature_directory,script_name,actual_directory):
	for file in os.listdir(os.path.join(feature_directory,"length_features")):
		df=pd.read_csv(os.path.join(feature_directory,"length_features",file),header=0)
		df=df.dropna()
		#print(df)
		df_list=[]
		names_id=df.ID.tolist()
		filename=file.replace("length_combined_", "")
		filename=filename.replace(".csv", "")
		for a in names_id:
			print(a)
			a_list=[]
			for b in names_id:
				command=os.popen("./profit -f {} -h {} {}".format(script_name, os.path.join(actual_directory,a+".pdb"), os.path.join(actual_directory,b+".pdb"))).readlines()
				if not os.path.exists(os.path.join(feature_directory,"length_RMSD_ByRes")):
					os.mkdir(os.path.join(feature_directory,"length_RMSD_ByRes"))
				By_res=open(os.path.join(feature_directory,"length_RMSD_ByRes",filename+ "_"+a+"_"+b+".profit_by_res"),'w+')
				counter=0
				list_write=[]
				for line in command:
					#print(line)
					columns=line.split("    ")
					if len(columns)< 2 and 'RMS:' in line:
						line=line.replace(" ","")
						line=line.replace("\n","")
						list_write.append(line[4:])
					if len(columns)>=3 and 'RMS:' in line:
						By_res.write(line)
						counter=1
					elif counter==1:
						By_res.write(" \n")
						counter=0
				if len(list_write)>7:
					del list_write[2]
					del list_write[3]
					del list_write[3]
					del list_write[4]
				if b==a:
					list_write=[0,0,0,0]

				list_write.append(a)
				list_write.append(b)
				if len(list_write) != 6:
					list_write=["Nan","Nan"]
				print(list_write)
				a_list=a_list+[list_write[1]]
			print("LOCAL_CA:", len(a_list)," ",len(names_id))
			df_list=df_list+[a_list]
		print(df_list)
		print(names_id)
		new_df=pd.DataFrame(df_list, columns=names_id, index=names_id)
		new_df.to_csv(os.path.join(feature_directory,"matrix_"+filename+".csv"), columns=names_id, index=names_id)



def delete_nans(matrix_directory):
	for file in os.listdir(os.path.join(matrix_directory)):
		df=pd.read_csv(os.path.join(matrix_directory,file),header=0, index_col=0)
		length=df.shape
		df2=df.replace(to_replace="Nan", value=np.nan, inplace=False)
		#df2=df2.dropna(axis=1, how='all', thresh=2)
		#df2=df2.dropna(axis=0, how='all', thresh=2)
		df2.to_csv(os.path.join(matrix_directory,file),header=True, index=False)



#sort_by_length(sys.argv[1],sys.argv[2])
#get_RMSDs(sys.argv[1],sys.argv[2],sys.argv[3])
#get_matrix(sys.argv[1],sys.argv[2],sys.argv[3])
delete_nans(sys.argv[1])