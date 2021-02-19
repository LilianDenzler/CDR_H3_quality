#!/usr/bin/env python
import qualiloop


import os


ID="test"
list_of_features=["loop_seq", "loop_length", "total_charge", "identity"]
file=str(os.path.join("/home/lilian/sync_project/CDRH3lib/tests","1F11_1.pdb"))
upload_dir="/var/tmp/"
print(myfunctions.merge(file, ID, list_of_features, upload_dir))


[['tip_pos', 'tip_res', 'protrusion', 'sequence', 'length',
       'total_charge', 'nr_charged', 'Happiness_mean', 'Nr_sad', 'identity',
       'similarity', 'template', 'target', 'Hydropathy', 'Hydropathy_diff',
       'Access', 'Relacc', 'Scacc', 'Screlacc', 'Access_avg', 'Relacc_avg',
       'Scacc_avg', 'Screlacc_avg', 'local_AA', 'global_AA', 'global_CA',
       'simlength', 'local_CA']]