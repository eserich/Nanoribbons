#!/usr/bin/python
## To create a nanoribbon of Graphene or hBN

import numpy as np
from math import acos 
from math import cos
from math import sqrt
from math import ceil

#Cosas que mejorar
#pasivation hBN
#que te diga el numero de atomos y el tamano de la celula
#xsf para vesta




############# INPUT

### Input structure file (txt file from vesta)
##Graphene
Gr_rect = np.array([[0.612500012, 0.707255504, 0.00000000],
                    [1.837500036, 1.414506708, 0.00000000],
                    [1.837500036, 2.829017717, 0.00000000],
                    [0.612500012, 3.536269111, 0.00000000]])

acell_Gr_rect = np.array([2.45,4.2435245514,0.0])

##hBN
hBN_rect = np.array([[0.621574998, 0.717737439, 0.00000000],
                     [1.864724994, 1.435470514, 0.00000000],
                     [1.864724994, 2.870941029, 0.00000000],
                     [0.621574998, 3.588678660, 0.00000000]])

acell_hBN_rect = np.array([2.4862999916,4.3064160347,0.0])

### Material
type_NR = raw_input('### Gr or hBN nanoribbon : ')


### How many cells

rows_in = raw_input('### Repetition of units cell in x : ')
rows    = int(rows_in.split()[0])


### Output file name
name = raw_input('### Ouput file name : ')


######### CREATE NEW STRUCTURE


def create_nanoribbon_structure(position_atom,acell,H_distance):
   global rows
   global name
   a = 0.5291772108
   pi= acos(-1)
   ##open output file
   w = open(name, 'w')
   
   print >>w, '%25s  %5i' %('## Armchair nanoribbon',rows)
   print >>w, '%25s  %10s' %('# Cartesian coordiantes','(Bohr)')
     
   H_reduction = (acell[0]/sqrt(3)-H_distance)*cos(pi/4)
   #H1
   print >>w, '%12.8f   %12.8f   %12.8f' %((position_atom[0,0]+H_reduction)/a , (position_atom[0,1]+H_reduction)/a , position_atom[0,2])  
   print >>w, '%12.8f   %12.8f   %12.8f' %(position_atom[1,0]/a , position_atom[1,1]/a , position_atom[0,2])
   print >>w, '%12.8f   %12.8f   %12.8f' %(position_atom[2,0]/a , position_atom[2,1]/a , position_atom[0,2]) 
   #H2
   print >>w, '%12.8f   %12.8f   %12.8f' %((position_atom[3,0]+H_reduction)/a , (position_atom[3,1]-H_reduction)/a , position_atom[0,2])

   i = 1			 
   if rows%2 == 0:
      while i <= (ceil(rows/2)-1):
         for k in range(4):
           print>>w, '%12.8f   %12.8f   %12.8f' %((position_atom[k,0]+i*acell[0])/a , position_atom[k,1]/a , position_atom[k,2])
         i+=1
      print >>w, '%12.8f   %12.8f   %12.8f' %((position_atom[0,0]+rows/2*acell[0])/a , (position_atom[0,1])/a , position_atom[0,2])
      print >>w, '%12.8f   %12.8f   %12.8f' %((position_atom[1,0]+rows/2*acell[0]-H_reduction)/a , (position_atom[1,1]-H_reduction)/a , position_atom[0,2])
      print >>w, '%12.8f   %12.8f   %12.8f' %((position_atom[2,0]+rows/2*acell[0]-H_reduction)/a , (position_atom[2,1]+H_reduction)/a , position_atom[0,2])
      print >>w, '%12.8f   %12.8f   %12.8f' %((position_atom[3,0]+rows/2*acell[0])/a , (position_atom[3,1])/a , position_atom[0,2])
   else:
       while i <= (ceil(rows/2)):
         for k in range(4):
           print>>w, '%12.8f   %12.8f   %12.8f' %((position_atom[k,0]+i*acell[0])/a , position_atom[k,1]/a , position_atom[k,2])
         i+=1
       print >>w, '%12.8f   %12.8f   %12.8f' %((position_atom[0,0]+(rows+1)/2*acell[0]-H_reduction)/a , (position_atom[0,1]+H_reduction)/a , position_atom[0,2])
       print >>w, '%12.8f   %12.8f   %12.8f' %((position_atom[3,0]+(rows+1)/2*acell[0]-H_reduction)/a , (position_atom[3,1]-H_reduction)/a , position_atom[0,2])

   w.close()



def write_vesta_file(position_atom,acell,H_distance):
    global rows
    global name

    name_xsf = name +'.xsf'
    w=open(name_xsf, 'w')

    print >>w, ' CRYSTAL'
    print >>w, ' PRIMVEC'
    print >>w, '%10.4f %10.4f %10.4f' %(acell[0]*round(rows/2)+8, 0.0000, 0.0000)
    print >>w, '%10.4f %10.4f %10.4f' %(0.0000 ,acell[1], 0.0000)
    print >>w, '%10.4f %10.4f %10.4f' %(0.0000, 0.0000, 14.0000)
    print >>w, ' PRIMCORD'
    print >>w, '%10i %10i' %(14,1)
    w.close() 







write_vesta_file(Gr_rect,acell_Gr_rect,1.07084)
			############ OUTPUT

if type_NR == 'Gr' :
   create_nanoribbon_structure(Gr_rect,acell_Gr_rect,1.07084)
elif type_NR == 'hBN' :
   create_nanoribbon_structure(hBN_rect,acell_hBN_rect,1.07084)
