#!/usr/bin/env python
#do chmod +x name_of_script before running

import pandas as pd

def get_charge(loop_seq_df):
	dic = {"D":-1, "K": -1,"R": 1,'E': 1, 'H':1}
	columns=["total_charge","nr_charged"]
	copy = False
	aminoacid=[]
	charge=0
	charged=0
	for index, row in loop_seq_df.iterrows():
		aminoacid=row["res"]
		if aminoacid in dic:
			charge+=dic[aminoacid]
			charged+=1
	write=[[charge,charged]]
	charge_df = pd.DataFrame(write, columns=columns, index=None)
	return (charge_df)