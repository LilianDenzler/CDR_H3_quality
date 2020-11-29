#!/usr/bin/env python
import sys
import os
#import numpy as np
import csv
import pandas as pd

def extract_seq(feature_directory,input_seq_directory):
    pass
def write_csv(feature_directory):
    pass


def get_loop_length (feature_directory,input_seq_directory):
    columns=["length","ID"]
    with open(os.path.join(feature_directory,"length" + ".csv"),"w") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
    for filename in os.listdir(input_seq_directory):
        name= filename.replace(".seq","")
        with open(os.path.join(input_seq_directory,filename)) as infile:
            copy = False
            count= 0
            for line in infile:
                if "H95" in line:
                    count+=1
                    copy = True
                    continue
                elif "H102" in line:
                    count+=1
                    copy = False
                    continue
                elif copy:
                    count+=1
        write=[str(count)]
        write.append(str(name))
        #print(write)
        with open(os.path.join(feature_directory,"length" + ".csv"),"a") as f:
            writer = csv.writer(f)
            writer.writerow(write)

def get_loop_charge(feature_directory,input_seq_directory):
    dic = {"D":-1, "K": -1,"R": 1,'E': 1, 'H':1}
    columns=["charge","nr_charged","ID"]
    with open(os.path.join(feature_directory,"charge" + ".csv"),"w") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
    for filename in os.listdir(input_seq_directory):
        name= filename.replace(".seq","")
        with open(os.path.join(input_seq_directory,filename)) as infile:
            copy = False
            aminoacid=[]
            charge=0
            nr_charged=0
            for line in infile:
                if "H95" in line:
                    copy = True
                if copy==True:
                    line2=line.split()
                    #print(line2)
                    aminoacid=line2[1]
                    #print ("AA:",aminoacid)
                    if aminoacid in dic:
                        charge+=dic[aminoacid]
                        nr_charged+=1
                        #print(name, charge, aminoacid, line2[0])
                if "H102" in line:
                    copy = False
        write=[charge,nr_charged,name]
        #print(write)
        with open(os.path.join(feature_directory,"charge" + ".csv"),"a") as f:
            writer = csv.writer(f)
            writer.writerow(write)

def happiness_score(feature_directory,actual_directory):
    '''Hydrophathy_index = {'A': 1.8, 'R': -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
                            "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2, "X": -0.5}#-0.5 is average
    columns=["Hydrophathy","ID"]
    with open(os.path.join(feature_directory,"hydropathy" + ".csv"),"w") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
    for filename in os.listdir(input_seq_directory):
        name= filename.replace(".seq","")
        hydro_value=0
        with open(os.path.join(input_seq_directory,filename)) as infile:
            copy = False
            for line in infile:
                if "H95" in line:
                    copy = True
                if copy==True:
                    line2=line.split()
                    aminoacid=line2[1]
                    pos=line2[0][1:]
                    hydro_value+=Hydrophathy_index[aminoacid]
                if "H102" in line:
                    copy = False
            write=[hydro_value,name]
            with open(os.path.join(feature_directory,"hydropathy_index_feature" + ".csv"),"a") as f:
                writer = csv.writer(f)
                writer.writerow(write)

        write=[hydro_value,name]
        with open(os.path.join(feature_directory,"hydropathy" + ".csv"),"a") as f:
            writer = csv.writer(f)
            writer.writerow(write)'''
    columns=["Happiness_mean","Nr_sad","ID"]
    with open(os.path.join(feature_directory,"happiness_score" + ".csv"),"w") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
    for filename in os.listdir(actual_directory):
        nr_sad=0
        name= filename.replace(".pdb","")
        command=os.popen("exposedhphob H95 H102 {}".format(os.path.join(actual_directory,filename))).readlines()
        for i in command:
            i=i.split()
            print (i)
            if "Mean:" in i[0]:
                mean=i[1]
                print(mean)
                continue
            if "Total:" in i[0]:
                pass
            elif float(i[2]) <0.5:
                nr_sad+=1
        write=[mean, nr_sad, name]
        with open(os.path.join(feature_directory,"happiness_score" + ".csv"),"a") as f:
            writer = csv.writer(f)
            writer.writerow(write)
    f.close()

def hydropathy(feature_directory, input_seq_directory, seq_id_file):
    '''Consensus values: Eisenberg, et al 'Faraday Symp.Chem.Soc'17(1982)109'''
    Hydrophathy_index = {'A': 00.250, 'R': -1.800, "N": -0.640, "D": -0.720, "C": 00.040, "Q": -0.690, "E": -0.620, "G": 00.160, "H": -0.400, "I": 00.730, "L": 00.530, "K": -1.100, "M": 00.260, "F": 00.610, "P": -0.070,
                            "S": -0.260, "T": -0.180, "W": 00.370, "Y": 00.020, "V": 00.540, "X": -0.5}#-0.5 is average
    columns_seq=["identity", "similarity", "template", "target", "ID"]
    data = pd.read_csv(sys.argv[3], names=columns_seq)
    #template=data.template.tolist()
    #target=data.target.tolist()


    columns=["Hydropathy","Hydropathy_diff","ID"]
    names=[]
    with open(os.path.join(feature_directory,"hydropathy" + ".csv"),"w") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
    for filename in os.listdir(input_seq_directory):
        name= filename.replace(".seq","")
        hydro_value=0
        with open(os.path.join(input_seq_directory,filename)) as infile:
            copy = False
            for line in infile:
                if "H95" in line:
                    copy = True
                if copy==True:
                    line2=line.split()
                    aminoacid=line2[1]
                    pos=line2[0][1:]
                    hydro_value+=Hydrophathy_index[aminoacid]
                    row=data.loc[data['ID'] == name]
                    if row.empty==True:
                        hydro_diff=None
                        continue
                    else:
                        template=row["template"].values[0]
                        target=row["target"].values[0]
                        for a,b in zip(template, target):
                            hydro_diff=abs(abs(Hydrophathy_index[b])-abs(Hydrophathy_index[a]))
                if "H102" in line:
                    copy = False
            write=[hydro_value,hydro_diff,name]
            with open(os.path.join(feature_directory,"hydropathy" + ".csv"),"a") as f:
                writer = csv.writer(f)
                writer.writerow(write)

        write=[hydro_value,hydro_diff, name]
        with open(os.path.join(feature_directory,"hydropathy" + ".csv"),"a") as f:
            writer = csv.writer(f)
            writer.writerow(write)   
         




def accessibility(feature_directory, actual_directory):
    '''Accessible_surface_area= {'A': 44.1, 'R': 152.9, "N": 80.8, "D": 76.3, "C": 56.4, "Q": 100.6, "E": 99.2, "G": 0, "H": 98.2, "I": 90.9, "L": 92.8, "K": 139.1, "M": 95.3, "F": 107.4, "P": 79.5,
                            "S": 57.5,"T": 73.4, "W": 143.4, "Y": 119.1, "V": 73, "X":89}#89 is rounded average
    columns=["Surface","ID"]
    with open(os.path.join(feature_directory,"accessible_surface" + ".csv"),"w") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
    for filename in os.listdir(input_seq_directory):
        name= filename.replace(".seq","")
        surface_values=[]
        with open(os.path.join(input_seq_directory,filename)) as infile:
            copy = False
            surface=0
            for line in infile:
                if "H95" in line:
                    copy = True
                if copy==True:
                    line2=line.split()
                    aminoacid=line2[1]
                    pos=line2[0][1:]
                    surface+=Accessible_surface_area[aminoacid]
                if "H102" in line:
                    copy = False
            write=[surface,name]         
            with open(os.path.join(feature_directory,"accessible_surface" + ".csv"),"a") as f:
                writer = csv.writer(f)
                writer.writerow(write)'''
    columns=["Access","Relacc","Scacc","Screlacc","Access_avg","Relacc_avg","Scacc_avg","Screlacc_avg","ID"]
    with open(os.path.join(feature_directory,"accessibility" + ".csv"),"w") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
    for filename in os.listdir(actual_directory):
        nr_sad=0
        name= filename.replace(".pdb","")
        command=os.popen("pdbsolv -r stdout {}".format(os.path.join(actual_directory,filename))).readlines()# RESIDUE  AA   ACCESS  RELACC  SCACC   SCRELACC
        copy=False
        results=False
        access=0
        relacc=0
        scacc=0
        screlacc=0
        counter=0
        for i in command:
            i=i.split()
            #print(i)
            if "END" in i:
                results=True
                continue
            if results==True:
                if "95"== i[2] and "H"==i[1]:
                    copy = True
                if copy==True:
                    counter+=1
                    access+=float(i[4])
                    relacc+=float(i[5])
                    scacc+=float(i[6])
                    screlacc+=float(i[7])
                if "102"==i[2] and "H"==i[1]:
                    copy = False
        if counter==0:
            access_avg=None
            relacc_avg=None
            scacc_avg=None
            screlacc_avg=None
            access=None
            relacc=None
            scacc=None
            screlacc=None
            write=[access,relacc,scacc,screlacc, access_avg,relacc_avg, scacc_avg, screlacc_avg, name]
            with open(os.path.join(feature_directory,"accessibility" + ".csv"),"a") as f:
                writer = csv.writer(f)
                writer.writerow(write)

        else:
            access_avg=access/counter
            relacc_avg=relacc/counter
            scacc_avg=scacc/counter
            screlacc_avg=screlacc/counter

            write=[access,relacc,scacc,screlacc, access_avg,relacc_avg, scacc_avg, screlacc_avg, name]
            with open(os.path.join(feature_directory,"accessibility" + ".csv"),"a") as f:
                writer = csv.writer(f)
                writer.writerow(write)
    f.close()



def zscale(feature_directory,input_seq_directory):
    pass
def VHSE_scale(feature_directory,input_seq_directory):
    pass
def seq_identity(feature_directory,input_seq_directory):
    pass
def bulged_non_bulged(feature_directory,input_seq_directory,model_directory,bioptools_directory):
    '''Previous studies identified a sequence motif (Arg or Lys at T2 and Asp at T6) that contributes to bulged torso formation in some but not all cases;
    these key residues were conserved in our bulged cluster, with 80% of bulged structures presenting Arg or Lys at T2, 73% presenting Asp at T6, and 65%
    retaining the complete T2/T6 sequence motif (S1 Fig) [9,10]. - Finn et al. PLOS One, 2016'''

    '''C-terminal loop residues form a pseudo dihedral angle,a101(100X–103 using Chothianumbering),
    of 39 ̊ and a pseudo bond angle,t101(100X–102using Chothia numbering), of 101 ̊

    They go on to use 3 standard deviations, so ca: pseudo dihedral a101=>6.5-70.7° and pseudo bond angle t101 =>81.81-117.81°
    t101, the Ca–Ca–Ca pseudo bond angle for the three C-terminal residues;and 2)a101,the Ca–Ca–Ca–Ca pseudo dihedral angle for
    the three C-terminal residues in the CDR H3loop and one adjacent residue inthe H chain framework

    -Weizner and Gray, 2016, Journal of Immunology'''

    '''for filename in os.listdir(input_seq_directory):
        file=open(os.path.join(input_seq_directory,filename))
        for line in file:
            line.split(" ")
            if line[1]=="86":'''
    pwd=os.getcwd()
    os.chdir(bioptools_directory)
    for filename in os.listdir(model_directory):
        name= filename.replace(".pdb.model","")
        angles=os.popen("pdbtorsions -c {}".format(os.path.join(model_directory,name+".pdb"+".model"))).readlines()
        for line in angles:
            columns=line.split()
            print (line)
            if "H100" in columns[0]:
                angle100x=columns[2]
                AAangle100x=columns[1]
                if float(angle100x) >=6.5 and float(angle100x) <= 70.7:
                    angle100x="k"
                else:
                    angle100x="e"
            if "H101" in columns[0]:
                angle101=columns[2]
                AAangle101=columns[1]
                if float(angle101) >=6.5 and float(angle101) <= 70.7:
                    angle101="k"
                else:
                    angle101="e"
            if "H102" in columns[0]:
                angle102=columns[2]
                AAangle102=columns[1]
                if float(angle102) >=6.5 and float(angle102) <= 70.7:
                    angle102="k"
                else:
                    angle102="e"
            if "H103" in columns[0]:
                angle103=columns[2]
                AAangle103=columns[1]
                if float(angle103) >=6.5 and float(angle103) <= 70.7:
                    angle103="k"
                else:
                    angle103="e"
                write=[angle100x, angle101, angle102, angle103, name]
                columnids=["angle100x","angle101","angle102","angle103","ID"]
                if os.path.exists(os.path.join(feature_directory,"kinked"+ ".csv")): #use feature_directory,"RMS_by_res_feature_csv_"+dics2[y] + ".csv") for separate files
                    with open(os.path.join(feature_directory,"kinked" + ".csv"), "a") as f:
                        writer = csv.writer(f)
                        writer.writerow(write)
                else:
                    with open(os.path.join(feature_directory,"kinked"+ ".csv"),"w") as f:
                        writer = csv.writer(f)
                        writer.writerow(columnids)
                        writer.writerow(write)


def nr_contacts(feature_directory,input_seq_directory):
    pass
def seq_ransomisation(feature_directory,input_seq_directory):
    pass
    #also information entropy
def BLOSUM_62():
    pass
def BLOSUM_50():
    pass

def similarity_score(feature_directory,log_directory):
    for filename in os.listdir(log_directory):
        name= filename.replace(".log","")
        with open(os.path.join(log_directory,filename)) as infile:
            for line in infile:
               #INFO: CDR-H3 (YEIR/YEWA) SeqID: 0.500 Similarity: 0.381
                if "INFO: CDR-H3" in line:
                    output=line.split(" ")
                    template_target=output[2]
                    template_target=template_target.split("/")
                    target=template_target[0]
                    template= template_target[1]
                    target=target.replace("(", "")
                    template=template.replace(")", "")

                    identity=output[4]
                    similarity=output[6]
                    similarity=similarity.replace('"', "")
                    similarity=similarity.replace("\n","")
                    write=[identity,similarity, template, target, name]
                    columnids=["identity","similarity", "template", "target", "ID"]
                    if os.path.exists(os.path.join(feature_directory,"seq_id"+ ".csv")): #use feature_directory,"RMS_by_res_feature_csv_"+dics2[y] + ".csv") for separate files
                        with open(os.path.join(feature_directory,"seq_id" + ".csv"), "a") as f:
                            writer = csv.writer(f)
                            writer.writerow(write)
                    else:
                        with open(os.path.join(feature_directory,"seq_id"+ ".csv"),"w") as f:
                            writer = csv.writer(f)
                            writer.writerow(columnids)
                            writer.writerow(write)


#also, how many AA of same category after one another
'''
Charged (side chains often form salt bridges):
• Arginine - Arg - R
• Lysine - Lys - K
• Aspartic acid - Asp - D
• Glutamic acid - Glu - E

Polar (form hydrogen bonds as proton donors or acceptors):
• Glutamine - Gln - Q
• Asparagine - Asn - N
• Histidine - His - H
• Serine - Ser - S
• Threonine - Thr - T
• Tyrosine - Tyr - Y
• Cysteine - Cys - C

Amphipathic (often found at the surface of proteins or lipid membranes, sometimes also classified as polar):
• Tryptophan - Trp - W
• Tyrosine - Tyr - Y
• Methionine - Met - M (may function as a ligand to metal ions)

Hydrophobic (normally buried inside the protein core):
• Alanine - Ala - A
• Isoleucine - Ile - I
• Leucine - Leu - L
• Methionine - Met - M
• Phenylalanine - Phe - F
• Valine - Val - V
• Proline - Pro - P
• Glycine - Gly - G
'''

def merge_csv(new_name, feature_directory):
    '''
    a = pd.read_csv(filea, header=0,sep=',')
    print(a)
    b = pd.read_csv(fileb, header=0,sep=',')
    print(b)
    #b = b.dropna(axis=1)
    merged = a.merge(b, on='ID')
    merged.to_csv(os.path.join(feature_directory,new_name+".csv"), index=False)'''
    file_list=[]
    for file in os.listdir(feature_directory):
        file_list.append(os.path.join(feature_directory,file))
    print(file_list)
    a = pd.read_csv(file_list[0], header=0,sep=',')
    for i in range (1,len(file_list)):
        bla=file_list[i]
        print(bla)
        b= pd.read_csv(str(bla), header=0,sep=',')
        merged = a.merge(b, on='ID')
        a=merged
    a.to_csv(os.path.join(feature_directory,new_name+".csv"), index=False)

#python3 feature_extraction.py ~/sync_project/Feature/ ~/sync_project/input_Abymod_seq/

if __name__=="__main__":
  #get_loop_length(sys.argv[1],sys.argv[2])
  #get_loop_charge(sys.argv[1],sys.argv[2])
  #bulged_non_bulged(sys.argv[1],sys.argv[2], sys.argv[3], sys.argv[4])
  #merge_csv(sys.argv[1], sys.argv[2],sys.argv[3], sys.argv[4])
  #similarity_score(sys.argv[1],sys.argv[2])
  #happiness_score(sys.argv[1],sys.argv[2])
  #hydropathy(sys.argv[1],sys.argv[2], sys.argv[3])
  #accessibility(sys.argv[1],sys.argv[2])
  feature_directory=sys.argv[1]
  new_name="all_combined_preliminary"
  merge_csv(new_name, feature_directory)
  