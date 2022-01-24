#!/usr/bin/python

# Obtain gap width from output file abinit

import numpy as np

DFT_EIG_in = raw_input(' ### insert output EIG file  :   ')
DFT_EIG    = open(DFT_EIG_in, 'r').readlines()


atoms_in  = raw_input(' ### Number and type of atoms  :   ')
atoms     = atoms_in.split()

ntype_atoms= len(atoms)/2




def obtain_number_occupated_bands(atoms):
    
    ntype_atoms = len(atoms)/2
    aux_vector  = np.zeros(ntype_atoms)
    for i in range(ntype_atoms):
        if   atoms[2*i+1] == 'H' : aux_vector[i]= 1*int(atoms[2*i])
        elif atoms[2*i+1] == 'B' : aux_vector[i]= 3*int(atoms[2*i])
        elif atoms[2*i+1] == 'C' : aux_vector[i]= 4*int(atoms[2*i])
        elif atoms[2*i+1] == 'N' : aux_vector[i]= 5*int(atoms[2*i])

    k=np.sum(aux_vector)
    return(k/2)

a=obtain_number_occupated_bands(atoms)


def obtain_gap(nEf_band, ndvik):
    global DFT_EIG
  
    print(float(DFT_EIG[int(2+(nEf_band-1)*(ndvik+2))+1].split()[1])-float(DFT_EIG[int(2+(nEf_band)*(ndvik+2))+1].split()[1]))

obtain_gap(a,51)
    





