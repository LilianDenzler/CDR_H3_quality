#!/usr/bin/env python
import sys
import os
import numpy as np
import csv
import pandas as pd
import itertools
import math
from biopandas.pdb import PandasPdb
from scipy.spatial.transform import Rotation


def get_line_len(ppdb):
	#P1=res H95, P2=res H102
	ppdb95=ppdb.loc[ppdb['residue_number'] == 95]
	ppdb102=ppdb.loc[ppdb['residue_number'] == 102]
	if ppdb95.empty or ppdb102.empty:
		return None
	P1x=float(ppdb95.x_coord)
	P1y=float(ppdb95.y_coord)
	P1z=float(ppdb95.z_coord)
	P2x=float(ppdb102.x_coord)
	P2y=float(ppdb102.y_coord)
	P2z=float(ppdb102.z_coord)
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

def get_PDB_coords(pdb_file):
	if pdb_file.endswith(".pdb")== False:
		print("this file is not in correct format (i.e. not PDB){}".format(filename))
	length_dic={}
	ppdb = PandasPdb().read_pdb(os.path.join(pdb_file))
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
		P=(Px,Py,Pz,R)
		#print(P)
		length_dic.update( {P : length} )
	max_value = max(length_dic.values())  # maximum value
	max_keys = [k for k, v in length_dic.items() if v == max_value]
	if len(max_keys)==1:
		max_keys=list(max_keys)
		max_keys=max_keys[0]
		#print(max_keys)
	else:
		max_keys==None

	P=[]
	R=[]
	for i in max_keys:
		P.append(i)
	R=list(P[3])
	P=P[:-1]
	return P,R,max_value

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
	R=(Rx,Ry,Rz)
	#print("frac",frac,"Rx",Rx,"Ry",Ry,"Rz",Rz, "P1x",P1x,"P2x",P2x)
	return length, R




def place_at_origin(pdb_file):
	if pdb_file.endswith(".pdb")== False:
		print("this file is not in correct format (i.e. not PDB){}".format(filename))
	ppdb = PandasPdb().read_pdb(os.path.join(pdb_file))
	ppdb=ppdb.df['ATOM']

	#get index of start-atom and stop-atom of loop
	index_95 = ppdb.loc[ppdb['residue_number'] == 95].index.tolist()[0]
	index_102 = ppdb.loc[ppdb['residue_number'] == 102].index.tolist()[-1]

	#get coordinates of full structure and of loop only
	full_structure=ppdb[["x_coord", "y_coord", "z_coord"]]
	full_structure=full_structure.to_numpy()
	loop_structure=ppdb.loc[ppdb['residue_number'] >= 95]
	loop_structure=loop_structure.loc[loop_structure['residue_number'] <= 102]
	loop_structure=loop_structure[["x_coord", "y_coord", "z_coord"]]
	
	#get base 95 point
	ppdb95=ppdb.loc[ppdb['residue_number'] == 95]
	ppdb95=ppdb95.loc[ppdb['chain_id'] == "H"]
	ppdb95=ppdb95.loc[ppdb['atom_name'] == "CA"]
	ppdb95=ppdb95[["x_coord","y_coord", "z_coord"]]
	ppdb95=ppdb95.to_numpy()

	#get base 102
	ppdb102=ppdb.loc[ppdb['residue_number'] == 102]
	ppdb102=ppdb102.loc[ppdb['chain_id'] == "H"]
	ppdb102=ppdb102.loc[ppdb['atom_name'] == "CA"]
	ppdb102=ppdb102[["x_coord","y_coord", "z_coord"]]
	ppdb102=ppdb102.to_numpy()

	#print(full_structure)
	#move entire strcuture so that resH95 is at origin
	for i in range(0, np.shape(full_structure)[0]):
		#print(i)
		#for x
		full_structure[i, 0]=full_structure[i, 0]-ppdb95[0,0]
		#for y
		full_structure[i, 1]=full_structure[i, 1]-ppdb95[0,1]
		#for z
		full_structure[i, 2]=full_structure[i, 2]-ppdb95[0,2]
	
	ppdb102[0,0]=ppdb102[0,0]-ppdb95[0,0]
	ppdb102[0,1]=ppdb102[0,1]-ppdb95[0,1]
	ppdb102[0,2]=ppdb102[0,2]-ppdb95[0,2]

	#align vectors from origin to res H102 and vector from origin along x axis
	rotation_structure, rmsd=Rotation.align_vectors(ppdb102,[[1,0,0]])

	#print(rotation_structure.as_matrix())

	rotated_full_structure=rotation_structure.apply(full_structure)
	rotated_full_df = pd.DataFrame(data=rotated_full_structure, columns=["x_coord", "y_coord", "z_coord"])
	#print(rotated_full_structure)

	print(rotated_full_df.iloc[500:1200])
	tip_point,close_point,max_length=get_PDB_coords(pdb_file)

	rotated_full_df=flap_loop(rotated_full_df, index_95,index_102)

	return (rotated_full_df, rotation_structure, ppdb102)


def flap_loop(rotated_full_df, index_95,index_102):
	# rotate the loop only for 1 degree around x axis
	r = Rotation.from_euler('x', 1, degrees=True)
	for i in range(index_95, index_102+1):
		a_loop_vector=rotated_full_df.iloc[[i]]
		new_loop_vector=r.apply(a_loop_vector)[0]
		rotated_full_df.loc[i] = new_loop_vector
	return (rotated_full_df)



if __name__ == '__main__':
  place_at_origin(sys.argv[1])
  