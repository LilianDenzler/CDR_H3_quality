#!/usr/bin/env python
import sys
import os
#import numpy as np
import csv
import pandas as pd

def check_loopdb(log_directory):
	list=[]
	for filename in os.listdir(log_directory):
		name= filename.replace(".log","")
		with open(os.path.join(log_directory,filename)) as infile:
			for line in infile:
				if "Running LoopDB" in line:
					list.append(name)
	print(list)
	return(list)
