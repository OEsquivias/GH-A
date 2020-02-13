#Get Puck program Python/numpy version
#Per començar: inputs el mateix que amb el getPuck original,
#per poder fer ràpides comprovacions

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
# sugar_traj.xyz=sugar_traj.xyz*10.0 #from nm to angstrom
sugar_xyz=sugar_traj.xyz*10.0 #from nm to angstrom

#Load the indexes from extract_sugars.tcl, put in array
# indexes=np.loadtxt("index.gt18",dtype=int)
# indexes=indexes-1 #Take into account it beginning in zero

#The number of atoms in the ring
atoms_ring=len(indexes)

#Put the coordinates in order
# def ordering(atoms_ring,sugar_xyz,frame,indexes):
	
# 	coord=np.zeros((atoms_ring,3),dtype=dt)
	
# 	for i in range(atoms_ring):
# 		coord[i,:]=sugar_xyz[frame,i+1,:]
# 	return coord



#Put the coordinates in order
# def ordering(atoms_ring,sugar_traj,frame,indexes):
	
# 	coord=np.zeros((atoms_ring,3),dtype=dt)
	
# 	for i in range(atoms_ring):
# 		coord[i,0]=sugar_traj.xyz[frame,indexes[i],0]
# 		coord[i,1]=sugar_traj.xyz[frame,indexes[i],1]
# 		coord[i,2]=sugar_traj.xyz[frame,indexes[i],2]
# #	coord=coord.astype('float64',casting='safe')
# 	return coord


#Geometric center function
def geocentr(atoms_ring,coord):

	Gcoord=np.zeros((atoms_ring,3),dtype=float)	
	trans=np.zeros((3),dtype=float)	

	for j in range(atoms_ring):
#		trans[0]=trans[0]+coord[j,0]/atoms_ring
#		trans[1]=trans[1]+coord[j,1]/atoms_ring
#		trans[2]=trans[2]+coord[j,2]/atoms_ring
		trans+=coord[j,:]/atoms_ring
	for j in range(atoms_ring):
#		Gcoord[j,0]=coord[j,0]-trans[0]
#		Gcoord[j,1]=coord[j,1]-trans[1]
#		Gcoord[j,2]=coord[j,2]-trans[2]
		Gcoord[j,:]=coord[j,:]-trans
	return Gcoord

#Lineal transformation of the vectors functions
def Xvector(atoms_ring,Gcoord):
	Xvect=np.zeros((3),dtype=float)
	for j in range(atoms_ring):
		Xvect[0]=Xvect[0]+Gcoord[j,0]*np.sin(2.0*np.pi*j/atoms_ring)
		Xvect[1]=Xvect[1]+Gcoord[j,1]*np.sin(2.0*np.pi*j/atoms_ring)
		Xvect[2]=Xvect[2]+Gcoord[j,2]*np.sin(2.0*np.pi*j/atoms_ring)

	modul=np.linalg.norm(Xvect)

	Xvect[0]=Xvect[0]/modul
	Xvect[1]=Xvect[1]/modul
	Xvect[2]=Xvect[2]/modul

	return Xvect
	
def Yvector(atoms_ring,Gcoord):
	Yvect=np.zeros((3),dtype=float)
	for j in range(atoms_ring):
                Yvect[0]=Yvect[0]+Gcoord[j,0]*np.cos(2.0*np.pi*j/atoms_ring)
                Yvect[1]=Yvect[1]+Gcoord[j,1]*np.cos(2.0*np.pi*j/atoms_ring)
                Yvect[2]=Yvect[2]+Gcoord[j,2]*np.cos(2.0*np.pi*j/atoms_ring)
	modulus=np.linalg.norm(Yvect)
	
	Yvect[0]=Yvect[0]/modulus
	Yvect[1]=Yvect[1]/modulus
	Yvect[2]=Yvect[2]/modulus

	return Yvect

def sumat1(atoms_ring,Zdispl):
	result=0
	for j in range(atoms_ring):
		result=result+Zdispl[j]*np.cos(2.0*np.pi/6.0*2.0*j)
	return result

def sumat2(atoms_ring,Zdispl):
	result=0
	for j in range(atoms_ring):
		result=result-Zdispl[j]*np.sin(2.0*np.pi/6.0*2.0*j)
	return result

def get_q2(atoms_ring,Zdispl):
	q2=np.sqrt(1.0/3.0)*np.sqrt(sumat1(atoms_ring,Zdispl)*sumat1(atoms_ring,Zdispl)+sumat2(atoms_ring,Zdispl)*sumat2(atoms_ring,Zdispl))
	return q2

def get_q3(atoms_ring,Zdispl):
	q3=np.sqrt(1.0/6.0)*(Zdispl[0]-Zdispl[1]+Zdispl[2]-Zdispl[3]+Zdispl[4]-Zdispl[5])
	return q3

def get_fi2(atoms_ring,Zdispl):
	fi2=np.arcsin(sumat2(atoms_ring,Zdispl)/np.sqrt(sumat1(atoms_ring,Zdispl)*sumat1(atoms_ring,Zdispl)+sumat2(atoms_ring,Zdispl)*sumat2(atoms_ring,Zdispl)))
	
	if sumat1(atoms_ring,Zdispl)<0:
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
for frame in range(sugar_traj.n_frames):  #to the number of frames in traj
#for frame in range(1):

#	print(sugar_traj.xyz[frame,:,:])

	#Put the coordinates in the correct order
	coord=ordering(atoms_ring,sugar_traj,frame,indexes)	
#	print(coord)

	#Get the centered coordinates
	Gcoord=geocentr(atoms_ring,coord)
#	print(Gcoord)

	#Get the new vectors Xvect,Yvect,Zvect - a partir d'aquí sembla que no està bé
	Xvect=Xvector(atoms_ring,Gcoord)
	Yvect=Yvector(atoms_ring,Gcoord)
	Zvect=np.zeros((atoms_ring,3))
	Zvect=np.cross(Xvect,Yvect)		
	
#	print(Xvect)
#	print(Yvect)
#	print(Zvect)

	#Calculate atomic displacements
	auxvect=np.zeros(3)
	Zdispl=np.zeros(atoms_ring,dtype=float)
	for atom in range(atoms_ring):
		auxvect[0]=Gcoord[atom,0]
		auxvect[1]=Gcoord[atom,1]
		auxvect[2]=Gcoord[atom,2]
		Zdispl[atom]=np.dot(auxvect,Zvect)

        #From this point,calculation of q2,q3 and fi2
#	print(frame+1)
#	print(get_q2(atoms_ring,Zdispl))	
#	print(get_fi2(atoms_ring,Zdispl)*360.0/2.0/np.pi)
#	print(get_q3(atoms_ring,Zdispl))
	q2=get_q2(atoms_ring,Zdispl)
	q3=get_q3(atoms_ring,Zdispl)
	fi2=get_fi2(atoms_ring,Zdispl)*360.0/2.0/np.pi

	#Calculation of Q and theta
	Q=get_q(q2,q3)
	theta=get_theta(Q,q3)
	print(Q,theta,fi2)	
