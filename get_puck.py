#Get Puck program Python/numpy version


#Canviar aquesta part perquè l'input sigui un fitxer de text amb aquesta informació
indexes=[8465,8467,8469,8471,8473,8480] #C1,C2,C3,C4,C5,O5
traj_file = "traj_1-81.dcd"
top_file = "gg_4glu_2.prmtop"

###

import numpy as np
import mdtraj as md

dt=np.dtype('f4')

#Load the trajectory from extract_sugars.tcl
# sugar_traj=md.load("gt18.xyz",top="gt18.pdb",)
sugar_traj=md.load(traj_file,top=top_file,atom_indices=indexes)
sugar_xyz=sugar_traj.xyz*10.0 #from nm to angstrom

#The number of atoms in the ring
global atoms_ring
atoms_ring=6


#Geometric center function
def geocentr(coord):
	centr=np.mean(coord, axis=0)
	gcoord=coord-centr
	return gcoord

#Lineal transformation of the vectors functions
def Xvector(Gcoord):
	Xvect=np.zeros((3),dtype=float)
	for j in range(atoms_ring):	
		Xvect[:]+=Gcoord[j,:]*np.sin(2.0*np.pi*j/atoms_ring)
	Xvect/=np.linalg.norm(Xvect)
	return Xvect

def Yvector(Gcoord):
	Yvect=np.zeros((3),dtype=float)
	for j in range(atoms_ring):	
 		Yvect[:]+=Gcoord[j,:]*np.cos(2.0*np.pi*j/atoms_ring)
	Yvect/=np.linalg.norm(Yvect)
	return Yvect

	#Calculate atomic displacements
def atom_dis(Gcoord):
	auxvect=np.zeros(3)
	Zdispl=np.zeros(atoms_ring,dtype=float)
	for atom in range(atoms_ring):
		auxvect=Gcoord[atom,:]
		Zdispl[atom]=np.dot(auxvect,Zvect)
	return Zdispl

def sumat1(Zdispl):
	result=0
	for j in range(atoms_ring):
		result+=Zdispl[j]*np.cos(2.0*np.pi*j/atoms_ring*2.0)
	return result
								
def sumat2(Zdispl):
	result=0
	for j in range(atoms_ring):
		result-=Zdispl[j]*np.sin(2.0*np.pi*j/atoms_ring*2.0)
	return result

def get_q2(Zdispl):
	q2=np.sqrt(1.0/3.0)*np.sqrt(sumat1(Zdispl)*sumat1(Zdispl)+sumat2(Zdispl)*sumat2(Zdispl))
	return q2

def get_q3(Zdispl):
	q3=np.sqrt(1.0/6.0)*(Zdispl[0]-Zdispl[1]+Zdispl[2]-Zdispl[3]+Zdispl[4]-Zdispl[5])
	return q3

def get_fi2(Zdispl):
	fi2=np.arcsin(sumat2(Zdispl)/np.sqrt(sumat1(Zdispl)*sumat1(Zdispl)+sumat2(Zdispl)*sumat2(Zdispl)))
	
	if sumat1(Zdispl)<0:
		fi2=np.pi-fi2
	elif fi2<0:
		fi2=2.0*np.pi+fi2
	
	return fi2
def get_q(q2,q3):
	q=np.sqrt((q2*q2)+(q3*q3))
	return q

def get_theta(q,q3):
	theta=np.arccos(q3/q)
	theta=theta*180.0/np.pi
	return theta




####################################################################################

#Main loop
output = np.zeros((sugar_traj.n_frames,3),dtype=float)
sugar_traj=md.load(traj_file,top=top_file,atom_indices=indexes)
for frame in range(sugar_traj.n_frames):  #to the number of frames in traj
	sugar_frame = sugar_traj.xyz[frame,:,:]*10.0


	#Get the centered coordinates
	Gcoord=geocentr(sugar_frame)

	#Get the new vectors Xvect,Yvect,Zvect
	Xvect=Xvector(Gcoord)
	Yvect=Yvector(Gcoord)
	Zvect=np.cross(Xvect,Yvect)		
	Zdispl=atom_dis(Gcoord)

    #From this point,calculation of q2,q3 and fi2
#	print(frame+1)
#	print(get_q2(Zdispl))	
#	print(get_fi2(Zdispl)*360.0/2.0/np.pi)
#	print(get_q3(Zdispl))
	q2=get_q2(Zdispl)
	q3=get_q3(Zdispl)
	fi2=get_fi2(Zdispl)*360.0/2.0/np.pi

	#Calculation of Q and theta
	Q=get_q(q2,q3)
	theta=get_theta(Q,q3)
	print(Q,theta,fi2)
	output[frame,:]=[Q,theta,fi2]
	# print(output[frame,:])
np.savetxt("output.out",output,fmt='%.8f',delimiter="\t")	
