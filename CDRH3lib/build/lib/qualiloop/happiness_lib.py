#!/usr/bin/env python
#do chmod +x name_of_script before running

import pandas as pd

def get_happy(file):
	columns=["Happiness_mean","Nr_sad"]
	nr_sad=0
	try:
		command=os.popen("exposedhphob H95 H102 {}".format(os.path.join(file))).readlines() 
		for i in command:
			i=i.split()
			if "Mean:" in i[0]:
				mean=i[1]
				#print(mean)
				continue
			if "Total:" in i[0]:
				pass
			elif float(i[2]) <0.5:
				nr_sad+=1
		write=[mean, nr_sad]
	except:
		print("happy: error ")

	happy_df = pd.DataFrame(write, columns=columns)
	return (happy_df)
