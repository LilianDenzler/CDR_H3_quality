#!/usr/bin/env python
import sys
import os
import numpy as np
import csv
#sys.argv[1]=input_seq_directory,
def loop_redundancy(input_seq_directory):
    list_seq={}
    doubles={}
    for file in os.listdir(input_seq_directory):
        with open(os.path.join(input_seq_directory,file)) as f:
            name=file.replace(".seq","")
            copy = False
            sequence=""
            for line in f:
                if "H95" in line:
                    copy = True
                if copy==True:
                    line2=line.split()
                    aminoacid=line2[1]
                    sequence+=str(aminoacid)
                if "H102" in line:
                    copy = False
            if sequence in list_seq:
                doubles.update({sequence:name})
            else:
                list_seq.update({sequence:name})
        f.close()
    for keys,values in doubles.items():
        doubles[keys]=[doubles[keys],list_seq[keys]]
    list_doubles1=[]
    list_doubles2=[]
    for keys, values in doubles.items():
        b=doubles[keys]
        a=b[0]
        list_doubles1.append(a)
    for keys, values in doubles.items():
        b=doubles[keys]
        a=b[1]
        list_doubles2.append(a)
    return (list_doubles1, list_doubles2)


def check_redundancy(to_be_checked_file, input_seq_directory):
    list_doubles1, list_doubles2= loop_redundancy(input_seq_directory)
    lines=[]
    with open(to_be_checked_file, 'r') as csvfile: 
        reader = csv.reader(csvfile)

        for row in reader:
            lines.append(row)
        for line in lines:
            for i in list_doubles1:
                if str(i) in line:
                    lines.remove(line)
                    #print("hey")

    with open('to_be_checked_file.csv', 'w') as writeFile:

        writer = csv.writer(writeFile)

        writer.writerows(lines)





if __name__ == '__main__':
    check_redundancy(sys.argv[1], sys.argv[2])

#just for checking:
'''for i in doubles.keys():
    for a in doubles[i]:
        with open(os.path.join(input_seq_directory,a+".seq")) as f:
            name=file.replace(".seq","")
            copy = False
            sequence=""
            for line in f:
                if "H95" in line:
                    copy = True
                if copy==True:
                    line2=line.split()
                    aminoacid=line2[1]
                    sequence+=str(aminoacid)
                if "H102" in line:
                    copy = False
            print(sequence)
        f.close()
'''
