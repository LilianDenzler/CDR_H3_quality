#!/usr/bin/env python
import qualiloop


import os


ID="test"
list_of_features=["loop_seq", "loop_length", "total_charge", "identity"]
file=str(os.path.join("/home/lilian/sync_project/CDRH3lib/tests","1F11_1.pdb"))
upload_dir="/var/tmp/"
print(myfunctions.merge(file, ID, list_of_features, upload_dir))