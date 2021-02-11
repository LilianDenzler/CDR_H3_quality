#!/usr/bin/env python
import sys
import os
import numpy as np
import csv
import pandas as pd
import itertools
import math
from biopandas.pdb import PandasPdb


def one_letter_code(residue):
  dic = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','XAA': 'X', 'UNK':'X'}
  #one_letter= atomium.structures.code()    #####doesnt work -> why?
  if len(residue) % 3 != 0:
    #print(residue)
    raise ValueError("error")
  one_letter= dic[residue]
  return (one_letter)

def parser(position,residue,chain_type):
  global string
  for i in zip(position, residue, chain_type):
    one_res=one_letter_code(residue)
    if chain_type=="H" or chain_type=="L":
      string+= '{}{}\t{}\n'.format(chain_type,position,one_res)
  return (string)

def get_line_len(ppdb):
  #print(ppdb)
  ppdb95=ppdb.loc[ppdb['residue_number'] == 95]
  ppdb105=ppdb.loc[ppdb['residue_number'] == 105]
  if ppdb95.empty or ppdb105.empty:
    return None
  P1x=float(ppdb95.x_coord)
  P1y=float(ppdb95.y_coord)
  P1z=float(ppdb95.z_coord)
  P2x=float(ppdb105.x_coord)
  P2y=float(ppdb105.y_coord)
  P2z=float(ppdb105.z_coord)
  #print(P1x,P1y,P1z,P2x,P2y,P2z)
  Lx=P2x-P1x
  Ly=P2y-P1y
  Lz=P2z-P1z
  line=[Lx,Ly,Lz]
  line_len=math.sqrt(np.dot(line,line))
  ux = Lx/line_len
  uy = Ly/line_len
  uz = Lz/line_len
  u=[ux,uy,uz]
  return line_len,P1x,P1y,P1z,P2x,P2y,P2z,u

def crossproduct(p1,p2):
  p3=[None,None,None]
  p3[0] = p1[1]*p2[2] - p1[2]*p2[1]
  p3[1] = p1[2]*p2[0] - p1[0]*p2[2]
  p3[2] = p1[0]*p2[1] - p1[1]*p2[0]
  return p3

def get_PDB_coords(filename, actual_directory):
  length_dic={}
  ppdb = PandasPdb().read_pdb(os.path.join(filename))
  ppdb=ppdb.df['ATOM']
  ppdb=ppdb.loc[ppdb['chain_id'] == "H"]
  ppdb=ppdb.loc[ppdb['atom_name'] == "CA"]
  if get_line_len(ppdb)==None:
    #print("getline")
    return None,None,None
  line_len,P1x,P1y,P1z,P2x,P2y,P2z,u=get_line_len(ppdb)
  ppdb=ppdb.loc[ppdb['residue_number'] > 95]
  ppdb=ppdb.loc[ppdb['residue_number'] < 105]
  #print(ppdb)
  for ind in ppdb.index: 
    Px=float(ppdb['x_coord'][ind])
    Py=float(ppdb['y_coord'][ind])
    Pz=float(ppdb['z_coord'][ind])
    P=[Px,Py,Pz]
    pos=ppdb["residue_number"][ind]
    length,R=calc(line_len,P,Px,Py,Pz,pos,P1x,P1y,P1z,P2x,P2y,P2z,u)
    length_dic.update( {pos : length} )
  max_value = max(length_dic.values())  # maximum value
  max_keys = [k for k, v in length_dic.items() if v == max_value]
  if len(max_keys)==1:
    max_keys=int(max_keys[0])
  else:
    max_keys==None
  ppdb_res=ppdb.loc[ppdb['residue_number'] == max_keys]
  res=ppdb_res.iloc[0]["residue_name"]
  res=one_letter_code(res)
  #print(length_dic)
  #print(max_keys)
  return max_keys,max_value,res

def calc(line_len,P,Px,Py,Pz,pos,P1x,P1y,P1z,P2x,P2y,P2z,u):
  frac=None
  Rx=None
  Ry=None
  Rz=None
  if line_len==0:
    length = math.sqrt((Px-P1x)*(Px-P1x) + (Py-P1y)*(Py-P1y) + (Pz-P1z)*(Pz-P1z))
    frac = 0;
    Rx = P1x
    Ry = P1y
    Rz = P1z
    print("same point")
    return length,R


  #Select Q as any point on line, we'll make it P1                      
  Qx = P1x
  Qy = P1y
  Qz = P1z
     
  #Calculate vector PQ                                               
  PQx = Qx - Px
  PQy = Qy - Py
  PQz = Qz - Pz
  PQ=[PQx,PQy,PQz]

  #Vector PR is the cross product of PQ and the unit vector along A (i.e. u)
  PR=crossproduct(PQ, u)
     
  #And the length of that vector is the length we want               
  length = math.sqrt(np.dot(PR,PR))
  #print("hey", length)
  #calculate where the closest point (R) on the line is to point P                                       

  #Find the projection of QP onto QP2                             
  QPx = Px - Qx
  QPy = Py - Qy
  QPz = Pz - Qz
  QP=[QPx,QPy,QPz]
      
  QP2x = P2x - Qx
  QP2y = P2y - Qy
  QP2z = P2z - Qz
  QP2=[QP2x,QP2y,QP2z]
      
  f = np.dot(QP, QP2) / math.sqrt(np.dot(QP2, QP2))
  if frac!=0:
    frac = f/line_len
  #Find point R: this is the fraction f of the unit vector along P1-P2 added onto Q 
  Rx = Qx + f * u[0]
  Ry = Qy + f * u[1]
  Rz = Qz + f * u[2]
  R=[Rx,Ry,Rz]
  #print("frac",frac,"Rx",Rx,"Ry",Ry,"Rz",Rz, "P1x",P1x,"P2x",P2x)
  return length, R
   
def calc_protrusion(filename):
  list_df=[]
  name=ID
  max_pos,max_length,res=get_PDB_coords(filename)
  line=[max_pos,res,max_length]
  list_df+=[line]
  df=pd.DataFrame(list_df,columns=["tip_pos","tip_res","protrusion"],index=None)
  return (df)


if __name__ == '__main__':
  df=calc_protrusion(sys.argv[1])
  df2=pd.read_csv(sys.argv[2],header=0)
  merged=df.merge(df2, on="ID")
  merged.to_csv(os.path.join(sys.argv[2],"new_merged.csv"),header=True,index=None)
