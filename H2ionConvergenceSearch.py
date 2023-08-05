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
    print("100th: ", EnT[100])
    print("200th: ", EnT[200])
    print("Number of energies: ", len(EnT))
    return EnT

#function that returns GE of a system
#interface function
def getGE(paraList):
    solve(paraList)
    return getEnergySpectrum(paraList)[0]

#this function outputs the number of digits 
#up to which the two floats are equal
#helper function
def compare_floats(float1, float2): 
    # Convert floats to strings
    str1 = str(float1)
    str2 = str(float2)
     
    # Remove negative signs
    if str1.startswith("-"):
        str1 = str1[1:]
    if str2.startswith("-"):
        str2 = str2[1:]

    # Remove the decimal point from the strings
    str1 = str1.replace(".", "")
    str2 = str2.replace(".", "")
    
    # Find the length of the shortest string
    min_length = min(len(str1), len(str2))
    
    # Find the index of the first difference
    first_difference = next((i for i, (digit1, digit2) in enumerate(zip(str1, str2)) if digit1 != digit2), min_length)
    
    return first_difference

#this function outputs parameterList, param in the parameterList has been
#updated to meet convergence criteria
#interface function
def findConvergence(paramName, parameters, incrementSteps, numDigits): 

    if not paramName in parameters:
        print("ERROR in function (findconvergence): paraName not in parameters")
        return parameters

    incrementStep = incrementSteps.get(paramName)

    #here checks convergence based on GE only, this can be improved later
    #TODO: improve this part for customisability
    currentGE = getGE(parameters)

    #if param is x_max -> here consider doubling x_max and N_x simultaneously
    #TODO: this obvisouly misses the "in-between" -> find a better way to find minimum basis set
    if paramName == "x_max": 
        while True: 
            #increment the x_max and N_x in parameterList simultaneously
            tempParameters = parameters
            tempParameters["x_max"] = parameters.get("x_max")*2 
            tempParameters["N_x"] = parameters.get("N_x")*2

            #compare new GE with previous GE
            tempGE = getGE(tempParameters)
            if compare_floats(currentGE, tempGE) >= numDigits: 
                parameters = tempParameters
                break
            
            #else update currentGE and parameters
            parameters = tempParameters
            currentGE = tempGE
            print("Current " + "x_max value: ", tempParameters.get("x_max"))
            print("Current " + "N_x value: ", tempParameters.get("N_x"))
    else: 
        while True: 
            #increment the param in parameterList
            tempParameters = parameters
            tempParameters[paramName] = parameters.get(paramName) + incrementSteps.get(paramName)

            #compare new GE with previous GE
            tempGE = getGE(tempParameters)
            if compare_floats(currentGE, tempGE) >= numDigits: 
                parameters = tempParameters
                break
            
            #else update currentGE and parameters
            parameters = tempParameters
            currentGE = tempGE
            print("Current " + paramName + " value: ", tempParameters.get(paramName))
    return (parameters)

#User input convergence criteria
convergenceDigits = 3
incrementSteps = {
    "N_x": 5,
    "N_y": 5,
    "k_x": 10,
    "k_y": 10
}


#helper function
def print_dict(dictionary):
    for key, value in dictionary.items():
        print(f'{key}: {value}')

#### TESTING EXECUTION ####

## here tests finding convergence for a single parameter ##
# print("starting parameters:")
# print_dict(parameters)
# getGE(parameters)
# print("finished running the starting parameters, now try finding convergence for N_x")
# findConvergence("N_x", parameters, incrementSteps, convergenceDigits)
# print("parameters after finding convergence for N_x:")
# print_dict(parameters)

## here tests finding convergence for all parameters ##
print("Finding convergence for all basis set parameters: ")
print("Starting parameters:")
print_dict(parameters)
print("Starting GE: ", getGE(parameters))

#QN: does the order by which I loop through these parameters matter?
basisParameters = ["x_max", "N_x", "N_y", "k_x", "k_y"]
for item in basisParameters:
    print("currently finding convergence for: " + item)
    findConvergence(item, parameters, incrementSteps, convergenceDigits)
    print("Finished finding convergence for: " + item)
    print("Current parameters:")
    print_dict(parameters)
    print("Curent GE: ", getGE(parameters))

print("Final parameters:")
print_dict(parameters)
print("Final GE: ", getGE(parameters))
