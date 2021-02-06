import sys
import os
import pandas as pd
#import save_RMS_lib
#import get_abYmod2_lib
#import seq_parser_lib
from pathlib import Path
import numpy as np
import epitopepredict as ep
from skimage.feature import hog
from skimage import data, exposure
from sklearn.decomposition import TruncatedSVD
import collections
from sklearn import manifold, datasets
from functools import partial

def get_input_seq(actual_directory, input_seq_directory):
    seq_parser_lib.run(actual_directory, input_seq_directory)

def get_models (redundant_file, input_seq_directory, model_directory):
    get_abYmod2_lib.run(redundant_file, input_seq_directory, model_directory)

def get_RMSDs (actual_directory, model_directory, feature_directory):
    df= save_RMS_lib.run(actual_directory, model_directory, feature_directory)
    print("get_RMSDs:",df.shape)
    return (df)

def get_loop_seq (input_seq_directory):
    columns=["sequence","ID"]
    data=[]
    for filename in os.listdir(input_seq_directory):
        name= filename.replace(".seq","")
        with open(os.path.join(input_seq_directory,filename)) as infile:
            copy = False
            seq=""
            for line in infile:
                if "H95" in line:
                    copy = True
                    seq+=(str((line.split())[1]))
                    continue
                elif "H102" in line:
                    seq+=(str((line.split())[1]))
                    copy = False
                    continue
                elif copy:
                    seq+=(str((line.split())[1]))
            write=[seq, name]
            data.append(write)
    df = pd.DataFrame(data, columns=columns)
    print("get_loop_seq:",df.shape)
    return (df)

def get_loop_length (input_seq_directory):
    columns=["length","ID"]
    data=[]
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
            write=[count,name]
            data.append(write)
    df = pd.DataFrame(data, columns=columns)
    print("get_loop_length:",df.shape)
    return (df)

def get_loop_charge(input_seq_directory):
    data=[]
    dic = {"D":-1, "K": -1,"R": 1,'E': 1, 'H':1}
    columns=["total_charge","nr_charged","ID"]
    for filename in os.listdir(input_seq_directory):
        name= filename.replace(".seq","")
        with open(os.path.join(input_seq_directory,filename)) as infile:
            copy = False
            aminoacid=[]
            charge=0
            charged=0
            for line in infile:
                if "H95" in line:
                    copy = True
                if copy==True:
                    line2=line.split()
                    aminoacid=line2[1]
                    if aminoacid in dic:
                        charge+=dic[aminoacid]
                        charged+=1
                if "H102" in line:
                    copy = False
            write=[charge,charged,name]
            data.append(write)
    df = pd.DataFrame(data, columns=columns)
    print("get_loop_charge:",df.shape)
    return (df)

def happiness_score(feature_directory,actual_directory):
    columns=["Happiness_mean","Nr_sad","ID"]
    data=[]
    for filename in os.listdir(actual_directory):
        nr_sad=0
        name= filename.replace(".pdb","")
        try:
            command=os.popen("exposedhphob H95 H102 {}".format(os.path.join(actual_directory,filename))).readlines()
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
            write=[mean, nr_sad, name]
            data.append(write)
        except:
            print("happy: error: ", name)

    df = pd.DataFrame(data, columns=columns)
    print("happiness_score:",df.shape)
    return (df)

def similarity_score(feature_directory,log_directory):
    columns=["identity","similarity", "template", "target", "ID"]
    data=[]
    for filename in os.listdir(log_directory):
        name= filename.replace(".log","")
        try:
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
                data.append(write)
        except:
            print("similarity: error: ", name)
    df = pd.DataFrame(data, columns=columns)
    print("similarity_score:",df.shape)
    return (df)
                    
def hydropathy(feature_directory, input_seq_directory,log_directory):
    columns=["Hydropathy","Hydropathy_diff","ID"]
    data=[]
    #Consensus values: Eisenberg, et al 'Faraday Symp.Chem.Soc'17(1982)109
    Hydrophathy_index = {'A': 00.250, 'R': -1.800, "N": -0.640, "D": -0.720, "C": 00.040, "Q": -0.690, "E": -0.620, "G": 00.160, "H": -0.400, "I": 00.730, "L": 00.530, "K": -1.100, "M": 00.260, "F": 00.610, "P": -0.070,
                            "S": -0.260, "T": -0.180, "W": 00.370, "Y": 00.020, "V": 00.540, "X": -0.5}#-0.5 is average
    df_sim=similarity_score(feature_directory, log_directory)
    for filename in os.listdir(input_seq_directory):
        name= filename.replace(".seq","")
        hydro_value=0
        with open(os.path.join(input_seq_directory,filename)) as infile:
            copy = False
            no_print=0
            for line in infile:
                if "H95" in line:
                    copy = True
                if copy==True:
                    line2=line.split()
                    aminoacid=line2[1]
                    pos=line2[0][1:]
                    hydro_value+=Hydrophathy_index[aminoacid]
                    row=df_sim.loc[df_sim['ID'] == name]
                    if row.empty==True:
                        hydro_diff=None
                        if no_print==0:
                            print(name)
                            no_print=1
                        continue
                    else:
                        template=row["template"].values[0]
                        target=row["target"].values[0]
                    for a,b in zip(template, target):
                        hydro_diff=abs(abs(Hydrophathy_index[b])-abs(Hydrophathy_index[a]))
                if "H102" in line:
                    copy = False
            write=[hydro_value,hydro_diff,name]
            data.append(write)
    df = pd.DataFrame(data, columns=columns)
    print("hydropathy:",df.shape)
    return (df)
                   
def accessibility(feature_directory, actual_directory):
    columns=["Access","Relacc","Scacc","Screlacc","Access_avg","Relacc_avg","Scacc_avg","Screlacc_avg","ID"]
    data=[]
    for filename in os.listdir(actual_directory):
        nr_sad=0
        name= filename.replace(".pdb","")
        try:
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
                print(i)
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
                data.append(write)

            else:
                access_avg=access/counter
                relacc_avg=relacc/counter
                scacc_avg=scacc/counter
                screlacc_avg=screlacc/counter

                write=[access,relacc,scacc,screlacc, access_avg,relacc_avg, scacc_avg, screlacc_avg, name]
                data.append(write)
        except:
            print("accessibility:error", name)
            
    df = pd.DataFrame(data, columns=columns)
    print("accessibility:",df.shape)
    return (df)


def merge(input_seq_directory, feature_directory, actual_directory, log_directory, model_directory):
    file_list=[]
    file_list.append(get_loop_seq(input_seq_directory))
    file_list.append(get_RMSDs(actual_directory, model_directory, feature_directory))
    file_list.append(get_loop_length(input_seq_directory))
    file_list.append(get_loop_charge(input_seq_directory))
    file_list.append(happiness_score(feature_directory,actual_directory))
    file_list.append(similarity_score(feature_directory,log_directory))
    file_list.append(hydropathy(feature_directory, input_seq_directory, log_directory))
    file_list.append(accessibility(feature_directory, actual_directory))
    file_list.append(RMSD_nom(actual_directory, model_directory, feature_directory))
    file_list.append(RMSD_binary(actual_directory, model_directory, feature_directory))
    a = file_list[0]
    for i in range (1,len(file_list)):
        b=file_list[i]
        merged = a.merge(b, on='ID')
        print("merged:", merged.shape)
        a=merged
    a.to_csv(os.path.join(feature_directory,"all_combined"+".csv"), index=False)
    print("final:", a.shape)
    all_combined = pd.read_csv(os.path.join(feature_directory,"all_combined"+".csv"))
    print("csv_final:", all_combined.shape)
    new_df=reduce(lambda df_left,df_right: pd.merge(df_left, df_right,left_index=True, right_index=True, how='outer'),file_list).fillna(nan_value)
    all_combined2 = pd.read_csv(os.path.join(feature_directory,"all_combined_newMethod"+".csv"))



def input1():
    actual_directory=input("Enter FULL path of directory containing all PDB files of actual structures:")
    actual_directory = Path(actual_directory)
    print(actual_directory)
    model_directory=input("Enter FULL path of directory containing all PDB files of structure models:")
    input_seq_directory=input("Enter FULL path of directory where you want to save your sequence files:")
    feature_directory=input("Enter FULL path of directory where you want to save your feature file:")
    log_directory=input("Enter FULL path of directory containing all log files of the modelling:")
    actual_directory = Path(actual_directory)
    model_directory = Path(model_directory)
    input_seq_directory = Path(input_seq_directory)
    feature_directory = Path(feature_directory)
    log_directory = Path(log_directory)
    get_input_seq(actual_directory, input_seq_directory)
    merge(input_seq_directory, feature_directory, actual_directory, log_directory, model_directory)

def input2():
    actual_directory=input("Enter FULL path of directory containing all PDB files of actual structures:")
    model_directory=input("Enter FULL path of directory where you want to save PDB files of structure models:")
    input_seq_directory=input("Enter FULL path of directory where you want to save your sequence files:")
    feature_directory=input("Enter FULL path of directory where you want to save your feature file:")
    log_directory=os.path.join(model_directory,"log_files")
    actual_directory = Path(actual_directory)
    model_directory = Path(model_directory)
    input_seq_directory = Path(input_seq_directory)
    feature_directory = Path(feature_directory)
    log_directory = Path(log_directory)
    get_input_seq(actual_directory, input_seq_directory)
    get_models (redundant_file, input_seq_directory, model_directory)
    merge(input_seq_directory, feature_directory, actual_directory, log_directory, model_directory)

def input3():
    actual_directory=input("Enter FULL path of directory containing all PDB files of actual structures:")
    model_directory=input("Enter FULL path of directory where you want to save PDB files of structure models:")
    input_seq_directory=input("Enter FULL path of directory containing your sequence files:")
    feature_directory=input("Enter FULL path of directory where you want to save your feature file:")
    log_directory=os.path.join(model_directory,"log_files")
    actual_directory = Path(actual_directory)
    model_directory = Path(model_directory)
    input_seq_directory = Path(input_seq_directory)
    feature_directory = Path(feature_directory)
    log_directory = Path(log_directory)
    get_models (redundant_file, input_seq_directory, model_directory)
    merge(input_seq_directory, feature_directory, actual_directory, log_directory, model_directory)



def input_type():
    correct=0
    what_input=input("What will you use to extract features of your models: \n 1. No modelling: PDB files of models & PDB files of actual structures (log files are present) \n 2. With modelling: PDB files of actual structures\n 3. With modelling: Input sequences of H and L chain \n TYPE: 1 / 2 / 3")
    if what_input=="1":
        input1()
    if what_input=="2":
        input2()
    if what_input=="3":
        input3()
    else:
        correct=1
        while correct==1:
            what_input=input("You did not input 1, 2, or 3. What will you use to extract features of your models:","\n","1. No modelling: PDB files of models & PDB files of actual structures",
        "\n","2. With modelling: PDB files of actual structures","\n","3. With modelling: Input sequences of H and L chain and PDB files of actual structures", "\n", "TYPE: 1 / 2 / 3")
            if what_input=="1":
                correct=0
            if what_input=="2":
                correct=0
            if what_input=="3":
                correct=0
            else:
                pass


def RMSD_binary(actual_directory, model_directory, feature_directory):
    columns=["local_AA_bin","local_CA_bin","global_AA_bin", "global_CA_bin", "ID"]
    data=[]
    threshold=int(input("Threshold value in ANgstrom (above-> not good model)"))
    RMSD_path=Path(input("enter FULL path of csv file containing RMSD values, press enter to calculate these"))
    if os.path.exists(RMSD_path) ==True:
        df=pd.read_csv(RMSD_path)
    else:
        input("AAAAH")
        df=get_RMSDs(actual_directory, model_directory, feature_directory)
    local_AA= df.local_AA.tolist()
    local_CA= df.local_CA.tolist()
    global_AA= df.global_AA.tolist()
    global_CA=df.global_CA.tolist()
    name=df.ID.tolist()
    for (a,b,c,d,n) in zip(local_AA, local_CA, global_AA, global_CA, name):
        a=float(a)
        b=float(b)
        c=float(c)
        d=float(d)
        if a>=threshold:
            a=1
        else:
            a=0
        if b>=threshold:
            b=1
        else:
            b=0
        if c>=threshold:
            c=1
        else:
            c=0
        if d>=threshold:
            d=1
        else:
            d=0
        write=[a,b,c,d,n]
        data.append(write)
    dfnew = pd.DataFrame(data, columns=columns)
    print("RMSD_binary:",dfnew.shape)
    return (dfnew)

def RMSD_nom(actual_directory, model_directory, feature_directory):
    columns=["local_AA_nom","local_CA_nom","global_AA_nom", "global_CA_nom", "ID"]
    data=[]
    RMSD_path=Path(input("enter FULL path of csv file containing RMSD values, press enter to calculate these"))
    if os.path.exists(RMSD_path) ==True:
        df=pd.read_csv(RMSD_path, header=0)
    else:
        input("AAAAH")
        df=get_RMSDs(actual_directory, model_directory, feature_directory)
    local_AA= df.local_AA.tolist()
    local_CA= df.local_CA.tolist()
    global_AA= df.global_AA.tolist()
    global_CA=df.global_CA.tolist()
    name=df.ID.tolist()
    threshlow=[0,2,4,6,8,10,12,14,16,18,20]
    threshup=[2,4,6,8,10,12,14,16,18,20,100]
    values=[1,2,3,4,5,6,7,8,9,10,11]
    print(len(threshlow), len(threshup), len(values))
    for (a,b,c,d,n) in zip(local_AA, local_CA, global_AA, global_CA, name):
        a=float(a)
        b=float(b)
        c=float(c)
        d=float(d)
        for (threshold_low,threshold_up,val) in list(zip(threshlow,threshup,values)):
            if type(a)==str:
                pass
            elif a>threshold_low and a<threshold_up:
                a=val
            if type(b)==str:
                pass
            elif b>threshold_low and b<threshold_up:
                b=val
            if type(c)==str:
                pass
            elif c>threshold_low and c<threshold_up:
                c=val
            if type(d)==str:
                pass
            elif d>threshold_low and d<threshold_up:
                d=val
        write=[int(a),int(b),int(c),int(d),n]
        data.append(write)
    dfnew = pd.DataFrame(data, columns=columns)
    print("RMSD_nom:",dfnew.shape)
    return (dfnew)

def one_hot(seq):
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L','M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    o = list(set(amino_acids) - set(seq))
    s = pd.DataFrame(list(seq))    
    x = pd.DataFrame(np.zeros((len(seq),len(o)),dtype=int),columns=o)    
    a = s[0].str.get_dummies(sep=',')
    a = a.join(x)
    a = a.sort_index(axis=1)
    #fd, hog_image = hog(a, visualize=True)
    #print(fd,hog_image)
    return (a)

def blosum62(seq):
    blosum = ep.blosum62
    seq=list(seq)
    x = pd.DataFrame([blosum[i] for i in seq]).reset_index(drop=True)
    #print(x)
    #la=input("hey")
    #x=list(x)
    #e = x.values.flatten()
    svd = TruncatedSVD(n_components=2)
    svd.fit(x)
    result = svd.transform(x)
    #print(result)
    blosum62_list=list(result)
    #fd, hog_image = hog(x, visualize=True)
    #print(fd,hog_image)
    #return result
    n_neighbors = 3
    n_components = 2
    LLE = partial(manifold.LocallyLinearEmbedding,n_neighbors, n_components, eigen_solver='auto')         #https://scikit-learn.org/stable/modules/manifold.html#manifold
    methods = collections.OrderedDict()
    methods['LLE'] = LLE(method='standard')
    methods['LTSA'] = LLE(method='ltsa')
    #methods['Hessian LLE'] = LLE(method='hessian')
    methods['Modified LLE'] = LLE(method='modified')
    methods['Isomap'] = manifold.Isomap(n_neighbors, n_components)
    methods['MDS'] = manifold.MDS(n_components, max_iter=100, n_init=1)
    methods['SE'] = manifold.SpectralEmbedding(n_components=n_components,n_neighbors=n_neighbors)        #https://scikit-learn.org/stable/auto_examples/manifold/plot_compare_methods.html
    methods['t-SNE'] = manifold.TSNE(n_components=n_components, init='pca',random_state=0)               #https://dmnfarrell.github.io/bioinformatics/mhclearning
    list_names=["LLE","LTSA", "Mod_LLE","Isomap", "MDS", "SE", "tSNE","blosum62"]
    dic_results=collections.OrderedDict()
    for i, (label, method) in enumerate(methods.items()):
        Y = method.fit_transform(x)
        dic_results[label]=Y
        #rint(label,Y)
        #a=input("HEYSS")
        #print(Y.shape)
    results_matrix=pd.DataFrame(zip(dic_results["LLE"],dic_results["LTSA"],dic_results["Modified LLE"],dic_results["Isomap"],dic_results["MDS"],dic_results["SE"],dic_results["t-SNE"],blosum62_list),columns=list_names)
    #print(results_matrix)
    return(results_matrix)


def critical_res(input_seq_directory):
    #Critical 6 residues (H94-H97, H101-H102)
    columns=["critical_seq","LLE","LTSA","Mod_LLE","Isomap","MDS","SE","tSNE","blosum62","ID"]
    data=[]
    for filename in os.listdir(input_seq_directory):
        name= filename.replace(".seq","")
        with open(os.path.join(input_seq_directory,filename)) as infile:
            copy = False
            for line in infile:
                if "H94" in line:
                    H94=str((line.split())[1])
                if "H95" in line:
                    H95=str((line.split())[1])
                if "H96" in line:
                    H96=str((line.split())[1])
                if "H97" in line:
                    H97=str((line.split())[1])
                if "H101" in line:
                    H101=str((line.split())[1])
                if "H102" in line:
                    H102=str((line.split())[1])
            seq=[H94,H95,H96,H97,H101,H102]
            blosum62_svd=blosum62(seq)
            LLE=blosum62_svd.LLE.tolist()
            LTSA=blosum62_svd.LTSA.tolist()
            Mod_LLE=blosum62_svd.Mod_LLE.tolist()
            Isomap=blosum62_svd.Isomap.tolist()
            MDS=blosum62_svd.MDS.tolist()
            SE=blosum62_svd.SE.tolist()
            tSNE=blosum62_svd.tSNE.tolist()
            blosum62_val=blosum62_svd.blosum62.tolist()
            seq_str=""
            seq_str=seq_str.join(seq)
            #print(seq_str)
            write=[seq_str,Mod_LLE,Isomap,blosum62_val]
            write_names=["critical_seq","critical_seq1", "critical_seq2","critical_seq3","critical_seq4", "critical_seq5", "critical_seq6",
            "Mod_LLE1","Mod_LLE2","Mod_LLE3","Mod_LLE4","Mod_LLE5","Mod_LLE6",
            "Mod_LLE1_2","Mod_LLE2_2","Mod_LLE3_2","Mod_LLE4_2","Mod_LLE5_2","Mod_LLE6_2",
            "Isomap1","Isomap2","Isomap3","Isomap4","Isomap5","Isomap6","Isomap1_2","Isomap2_2","Isomap3_2","Isomap4_2","Isomap5_2","Isomap6_2",
            "blosum62_val1","blosum62_val2","blosum62_val3","blosum62_val4","blosum62_val5","blosum62_val6",
            "blosum62_val1_2","blosum62_val2_2","blosum62_val3_2","blosum62_val4_2","blosum62_val5_2","blosum62_val6_2","ID"]
            new_write=[]
            new_write+=[seq_str]
            for i in write:
                print (i)
                for a in i:
                    print(a)
                    for b in a:
                        new_write.append(b)
            new_write.append(name)

            data.append(new_write)
            #bla=input(data)
    df = pd.DataFrame(data, columns=write_names)
    print("get_loop_seq:",df.shape)
    print(df)
    return (df)

    #python3 myfunctions.py ~/sync_project/input_Abymod_seq/

    #6x20 Integer (1-hot)6x20 
    #Real (BLOSUM80)6x20 
    #Real (BLOSUM62)6x20 
    #Real (BLOSUM45)6x20 
    #Real (PAM250/MDM)6x5 
    #Real (Li-Koehl PCA5)6x4 
    #Real (Li-Koehl PCA4)6x3 
    #Real (Li-Koehl PCA3)6x3 
    #Real (El-Maarty 3-parameter physical encoding)6x4 
    #Real (Abhinandan 4-parameter physical encoding)

def similarity_length(feature_file, feature_directory):
    df=pd.read_csv(feature_file, header=0)
    df1=df.copy()
    length=df.length.tolist()
    similarity=df.similarity.tolist()
    list1=[]
    for (a,b) in zip(length, similarity):
        list1=list1+[b/a]
    print(list1)
    df["simlength"]=list1
    column_names=list(df1.columns) 
    merged = df.merge(df1, on=column_names)
    print("merged:", merged.shape)
    merged.to_csv(os.path.join(feature_directory,"new"+".csv"), index=False)
    print("final:", merged.shape)

    #Sum and mean of (Hydrophobicity x residue solvent accessibility)

if __name__ == '__main__':
    #get_models(sys.argv[6], ssys.argv[1], sys.argv[5])
    #merge(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],sys.arvg[5])

    #input arguments: 
    #input_seq_directory, feature_directory, actual_directory, log_directory,model_directory,redundant_file

    #similarity_length(sys.argv[1], sys.argv[2])
    #input_type()
    #blosum62(["A","G","H", "L", "M"])
    feature_directory=sys.argv[3]
    file2=os.path.join(sys.argv[2])
    file2=pd.read_csv(file2,header=0)
    bla=input(file2)
    df1=critical_res(sys.argv[1])
    
    merged = df1.merge(file2, on='ID')
    print("merged:", merged.shape)
    a=merged
    a.to_csv(os.path.join(feature_directory,"features+critical"+".csv"), index=False)
    print("final:", a.shape)
