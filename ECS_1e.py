
#   cd /users/amo/walterni/h2_vibwavepacket/  && python H2wavepacket.py 0

# 0 = create Tasks
# 1 = see Outoput of B_DAM_ECS_1e
# 2 = see Outoput of B_DAM_ECS_2e
# 2 = see Outoput of Dipole_2e


import math

from math import e
import numpy as np
from scipy.io import FortranFile
import sys 
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline
from scipy.io import _fortran
import argparse
import sys
#import pandas as pd
import re
import os.path
from matplotlib import animation
import seaborn as sns
import pandas as pd
import numpy as np

#from UnitConversion import *

sys.path.insert(1, '/users/amo/walterni/MyGeneralScripts/')
from UnitConversion import *

Task=int(sys.argv[1])  # first: Number os task    #second:'task' or 'out' 


if Task==-1:
    sys.exit()



atom1='H'
atom2='H'


base='BaseC_R0p400' # info on: R, m_min, m_max, Pot-Type, Pot-Sym, state-selection, charge atoms
R_int=0.4
x_max=10
k_x,k_y=2,2
m_max=3




# info for 1e-output
hm_range=5
hm_points=50
nn=6

name=''.join((base,'_stateNr_',str(nn)))

E_int=1/R_int



##### TASKS
if Task==0:
    print('')
    print(' '.join(('d1 && RUN_B_DAM_ECS_1e.csh',atom1,atom2,base)))  
    print('')
    sys.exit()



### output paths
path0='/users/amo/walterni/AMO_TOOLS_USER/ELECTR_STR/B_DAM_ECS'  # = d1
folder_1e_out=''.join((path0,'/B_DAM_ECS_1e/out/',atom1,'_',atom2,'/',atom1,'_',atom2,'_',base))



#######################################################################################################
#######################################################################################################
#######################################################################################################
# #####################################################################################################






xx=np.linspace(1,x_max, 1000)
yy=np.linspace(-1,1, 1000)


# See OUTOUT OF RUN B_DAM_ECS_1e 


fig, axs= plt.subplots(2,m_max+1)

for symind in (0,1):
    for m in range(m_max+1):


        if symind==0:
            sym='g'
        if symind==1:
            sym='u'
        QZ_m=m



        pout1=''.join((folder_1e_out,'/m__',str(QZ_m),'_',sym,'.out1e')) 
        f_out1= FortranFile(pout1, 'r')
        n_calc_state=f_out1.read_ints(np.int32)[0]
        print('n_calc_state:', n_calc_state)
        EnT=[0]*n_calc_state
        CxT=[0]*n_calc_state
        CyT=[0]*n_calc_state
        for n in range(n_calc_state):
            EnT[n]=f_out1.read_reals(float)[0] #,why is this a complex value? second value small
            CxT[n]=f_out1.read_reals(float)    # Length B_N_X -1  (last 0ne is zero) 
            CyT[n]=f_out1.read_reals(float)     # Length B_N_y -0   
            if n==0:
                N_x,N_y=len(CxT[n])+1,len(CyT[n])
                print('N_x,N_y:', N_x,N_y)
        f_out1.close

        knots_x=np.concatenate(([1]*(k_x-1),np.linspace(1,x_max,N_x-k_x+2), [x_max]*(k_x-1)))
        knots_y=np.concatenate(([-1]*(k_y-1),np.linspace(-1,1,2*N_y-k_y+2), [1]*(k_y-1)))

        def bsp_x(x_coord,c):    #input:b-spline parameters of length B_N-2
            c0=np.concatenate((c,[0]))
            bsp = BSpline(knots_x, c0, k_x-1)
            return np.array(bsp(x_coord))  #output: y-values  
    
        def bsp_y(y_coord,c,gu):    #input:b-spline parameters of length B_N-2
            if gu=='g':
                c0=np.concatenate((c,c[::-1]))  # was guessed 
            if gu=='u':
                c0=np.concatenate((-c,c[::-1]))  # was guessed 
            bsp = BSpline(knots_y, c0, k_y-1)
            return np.array(bsp(y_coord))  #output: y-values 

        # fig,axs= plt.subplots(1,2)
        # fx=bsp_x(xx,CxT[nn])    
        # fx=fx*sign(fx[0])
        # axs[0].plot(xx,fx)
        # axs[0].set_xscale('log')
        # fy=bsp_y(yy,CyT[nn],sym)    
        # fy=fy*sign(fy[0])
        # axs[1].plot(yy,fy)


        def elec(z_,r_):
            xi,yi=cart2sphr(z_,r_,R_int)
            phi=bsp_x(xi,CxT[nn])*bsp_y(yi,CyT[nn],sym)
            return phi



        hm_zz = np.linspace(-hm_range*R_int, hm_range*R_int, hm_points)
        hm_rr = np.linspace(-hm_range*R_int, hm_range*R_int, hm_points)
        M=np.zeros((len(hm_zz),len(hm_rr)))
        for r in range(len(hm_rr)):
            for z in range(len(hm_zz)):
                M[r,z]= elec(hm_zz[z],hm_rr[r])**2

        if sym=='g':        
            sns.heatmap(M, cmap='Blues',xticklabels=False,yticklabels=False,cbar=False,ax=axs[symind,m])
        if sym=='u':       
            sns.heatmap(M, cmap='Greens',xticklabels=False,yticklabels=False,cbar=False,ax=axs[symind,m])

if Task==1:
    plt.show()


if Task==5:
    plt.savefig(''.join(('results/ECS_1e_',name,'.pdf')))
    print(name,'exported')
    plt.show()

