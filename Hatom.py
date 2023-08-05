import math
import scipy
from math import e
import numpy as np
from scipy.io import FortranFile
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline
from scipy.io import _fortran
import argparse
import re
import os.path
import os
from matplotlib import animation
import seaborn as sns
import pandas as pd
import numpy as np
from sys import argv #sys.argv is to input additional arguments to the python script when running it using command line in terminal 
import fortranformat as ff
#sys.path.insert(1, '/users/stud/linzhonglin/linzhonglin/MyGeneralScripts/')
from UnitConversion import *
from LaserExcitation import *
from File_reader import *
import shutil
#from termcolor import colored
#from matplotlib import animation
from scipy.integrate import solve_ivp
from Templates import *

Task=sys.argv[1] 
if Task=='z':    #use this to test output from other scripts 
    sys.exit()


# parameters
r_max = 100
k_r = 10
N_r = 600
max_l = 5
NameOfBase = 'BSplines'

#########################################################################################
######################### create files ##################################################
#########################################################################################

if Task=='ecs' and sys.argv[2].find('c')!=-1: 
    path_input_base='/users/stud/linzhonglin/AMO_TOOLS_USER/ELECTR_STR/HSPH/input/'  # Base
    #path_input_dip='/users/amo/walterni/AMO_TOOLS_USER/ELECTR_STR/B_DAM_ECS/DIPOLE_2e/input/H_H/'  #Dip
    filename=''.join((NameOfBase,str(r_max),'.hs'))
    pathname=''.join((path_input_base,filename))
    content_exp=createHinput(r_max, k_r, N_r, max_l)
    f_exp1 = open(pathname,'w')
    f_exp1.write(content_exp)
    f_exp1.close()
    print('created:',pathname)
    print("finished task c")

     
        
     

#########################################################################################
######################### run code files ##################################################
#########################################################################################

if Task=='ecs' and sys.argv[2].find('r')!=-1: 

    pathcode='/users/stud/linzhonglin/AMO_TOOLS_USER/ELECTR_STR/HSPH/'

    filename=''.join((NameOfBase,str(r_max)))
    ecs1e_task=' '.join((pathcode+'RUN_HSPH_DIAG.csh',filename))
    print(ecs1e_task)
    os.system(ecs1e_task)
    print("finished task r")


if Task=='ecs' and sys.argv[2].find('s')!=-1: 
    path_out_1e='/users/stud/linzhonglin/AMO_TOOLS_USER/ELECTR_STR/HSPH/out/'
    outputfilename=''.join((NameOfBase,str(r_max)))
    path_file=path_out_1e+outputfilename+'.ec'

    print("showing file:"+path_file)
    fout=FortranFile(path_file,'r')

    n_calc_state=fout.read_ints(np.int32)
    print('n_calc_state:', n_calc_state)
    m_calc_state=fout.read_ints(np.int32)
    print('n_calc_state:', n_calc_state)



    EnT=fout.read_reals(float)
    print('Energies:',EnT,len(EnT))
    n_calc_state=fout.read_reals(float)
    
    
    CrT=fout.read_reals(float)
    print('C-coefficients:',CrT,len(CrT))

    CrT=fout.read_reals(float)
    #print('C-coefficients:',CrT,len(CrT))
   


    # fig, axs= plt.subplots(2)
    axs[0].plot([0]*len(EnT),EnT,'o')
    plt.show()




   
    sys.exit()
    CxT[n]=fout.read_reals(float)    # Length B_N_X -1  (last 0ne is zero) 



    sys.exit()
    EnT=[0]*n_calc_state
    CxT=[0]*n_calc_state
    CyT=[0]*n_calc_state
    for n in range(n_calc_state):
        EnT[n]=fout.read_reals(float)[0] #,why is this a complex value? second value small

        sys.exit()
        CxT[n]=fout.read_reals(float)    # Length B_N_X -1  (last 0ne is zero) 
        CyT[n]=fout.read_reals(float)     # Length B_N_y -0   
        
        #since csr stps can be run separatly,
        #this if block is to ensure that when basis parameters are changed, 
        #r step should be run again
        if n==0:
            N_x_f,N_y_f=len(CxT[n])+1,len(CyT[n])
            print('N_x_f,N_y_f(from fortran output):', N_x_f,N_y_f)
            if (N_x_f != N_x) or (N_y_f != N_y):
                print('ERROR!!!! MISMATCH IN BASE')
                sys.exit()

    fout.close()

    print("000th: ", EnT[0])
    print("100th: ", EnT[100])
    print("200th: ", EnT[200])
    print(len(CxT[0]),len(CyT[0]))

    print("Number of energies: ", len(EnT))




    if sys.argv[2].find('w')!=-1: 
        knots_x=np.concatenate(([1]*(k_x-1),np.linspace(1,x_max,N_x-k_x+2), [x_max]*(k_x-1)))
        knots_y=np.concatenate(([-1]*(k_y-1),np.linspace(-1,1,2*N_y-k_y+2), [1]*(k_y-1)))
        print(knots_x)
        print(knots_y)

        def bsp_x(x_coord,c):    #input:b-spline parameters of length B_N-2
                c0=np.concatenate(([0],c))
                bsp = BSpline(knots_x, c0, k_x-1)
                return np.array(bsp(x_coord))  #output: y-values  
        
        def bsp_y(y_coord,c,gu):    #input:b-spline parameters of length B_N-2
           
            if gu=='g':
                c0=np.concatenate((c,c[::-1]))  # was guessed 
            if gu=='u':
                c0=np.concatenate((-c,c[::-1]))  # was guessed 


            bsp = BSpline(knots_y, c, k_y-1)
            return np.array(bsp(y_coord))  #output: y-values 

        def elec(z_,r_):
                    xi,yi=cart2sphr(z_,r_,R_int)
                    phi=bsp_x(xi,CxT[nn])*bsp_y(yi,CyT[nn],parity_show)
                    return phi
        

        fig, axs= plt.subplots(2,m_max+1)

        hm_range=5 #plot range
        hm_points=70 #plot resolution

        hm_zz = np.linspace(-hm_range*R_int, hm_range*R_int, hm_points)
        hm_rr = np.linspace(-hm_range*R_int, hm_range*R_int, hm_points)

        M=np.zeros((len(hm_zz),len(hm_rr)))
        for r in range(len(hm_rr)):
            for z in range(len(hm_zz)):
                M[r,z]= elec(hm_zz[z],hm_rr[r])**2

        if parity_show=='g':        
            sns.heatmap(M, cmap='Blues',xticklabels=False,yticklabels=False,cbar=False,ax=axs[0,0])
        if parity_show=='u':       
            sns.heatmap(M, cmap='Greens',xticklabels=False,yticklabels=False,cbar=False,ax=axs[0,0])

        plt.show()




