#!/usr/bin/env python

import os
import sys

for filename in os.listdir(sys.argv[1]):
	structure=(sys.argv[1]+"/"+filename)
  subprocess = subprocess.Popen("profit", "-f", "RMS_calc.txt", structure, shell=True, stdout=subprocess.PIPE)
  return_code = subprocess.poll()
  if return_code is not None:
    subprocess_return = subprocess.stdout.read()
