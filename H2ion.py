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
from UnitConversion import *

#these are starting parameters
parameters = {
    "R_int": 10,
    "N_dig": 3,
    "m_min": 1,
    "m_max": 1,
    "NameOfBase": 'Bspline', #if NameOfBase let empty, it will be called according to the number of AO
    "x_max": 1000,
    "N_x": 30,
    "N_y": 30,
    "k_x": 10,
    "k_y": 10
}


#solving for energy and wavefunction
parity_show='u'
m_show=1
nn=100

#sovlve the system with basis as given by input parameters through FORRTAN code
#helper function
def solve(parameterList):
    R_int=parameterList.get("R_int")
    N_dig=parameterList.get("N_dig")
    m_min=parameterList.get("m_min")
    m_max=parameterList.get("m_max")
    NameOfBase=parameterList.get("NameOfBase")
    x_max=parameterList.get("x_max")
    N_x=parameterList.get("N_x")
    N_y=parameterList.get("N_y")
    k_x=parameterList.get("k_x")
    k_y=parameterList.get("k_y")

    # create input file for fortan code 
    path_input_base='/users/stud/linzhonglin/AMO_TOOLS_USER/ELECTR_STR/B_DAM_ECS/B_DAM_ECS_1e/input/H_H/'  # Base
    filename=''.join(('H_H_',NameOfBase,'_R',p_str(R_int,N_dig),'.bde1'))
    pathname=''.join((path_input_base,filename))
    content_exp=bde1(float(R_int),m_min,m_max,x_max,N_x,N_y,k_x,k_y)
    f_exp1 = open(pathname,'w')
    f_exp1.write(content_exp)
    f_exp1.close()
    # print('created FORTRAN input file:',pathname)
    
    #run the fortan code with the created input file
    pathcode='/users/stud/linzhonglin/AMO_TOOLS_USER/ELECTR_STR/B_DAM_ECS/'
    inputfile=''.join((NameOfBase,'_R',p_str(R_int,N_dig)))
    ecs1e_task=' '.join((pathcode+'B_DAM_ECS_1e/RUN_B_DAM_ECS_1e.csh','H','H',inputfile))
    # print(ecs1e_task)
    os.system(ecs1e_task) #this runs the FORTRAN code through terminal (probably)

#reading the energy spectrum from FORTRAN output file (for the given set of parameters)
#interface function
def getEnergySpectrum(parameterList):
    R_int=parameterList.get("R_int")
    N_dig=parameterList.get("N_dig")
    m_min=parameterList.get("m_min")
    m_max=parameterList.get("m_max")
    NameOfBase=parameterList.get("NameOfBase")
    x_max=parameterList.get("x_max")
    N_x=parameterList.get("N_x")
    N_y=parameterList.get("N_y")
    k_x=parameterList.get("k_x")
    k_y=parameterList.get("k_y")

    solve(parameterList)

    path_out_1e='/users/stud/linzhonglin/AMO_TOOLS_USER/ELECTR_STR/B_DAM_ECS/B_DAM_ECS_1e/out/H_H/H_H_'
    path_file=path_out_1e+NameOfBase+"_R"+p_str(R_int,N_dig)+"/m__"+str(m_show)+"_"+parity_show+".out1e"
    path_file=path_out_1e+NameOfBase+"_R"+p_str(R_int,N_dig)+"/m__"+str(m_show)+".out1e"
    # print("showing file:"+path_file)
    fout=FortranFile(path_file,'r')

    n_calc_state=fout.read_ints(np.int32)[0]
    # print('n_calc_state:', n_calc_state)

    EnT=[0]*n_calc_state
    CxT=[0]*n_calc_state
    CyT=[0]*n_calc_state
    for n in range(n_calc_state):
        EnT[n]=fout.read_reals(float)[0] #,why is this a complex value? second value small
        CxT[n]=fout.read_reals(float)    # Length B_N_X -1  (last 0ne is zero) 
        CyT[n]=fout.read_reals(float)     # Length B_N_y -0   
    
    fout.close()
    print("000th: ", EnT[0])
    print("Cx-coefficients:", CxT[0])
    print("Number of Cx-coefficients:", len(CxT[0]))
    print("Cy-coefficients:", CxT[0])
    print("Number of Cy-coefficients:", len(CyT[0]))
    # print("100th: ", EnT[100])
    # print("200th: ", EnT[200])
    print("Number of energies: ", len(EnT))


    def bsp_x(x_coord,c):    #input:b-spline parameters of length B_N-2
                c0=np.concatenate((c,[0]))
                bsp = BSpline(knots_x, c0, k_x-1)
                return np.array(bsp(x_coord))  #output: y-values  
    
    
    return EnT

#function that returns GE of a system
#interface function
def getGE(paraList):
    solve(paraList)
    return getEnergySpectrum(paraList)[0]


#### TESTING EXECUTION ####
getEnergySpectrum(parameters)



