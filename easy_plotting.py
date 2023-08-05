#%%
import plotting
import importlib
importlib.reload(plotting)
import parameters
importlib.reload(parameters)
import matplotlib.pyplot as plt
import ECS12
import run_multiple_parameter

system_basis_parameters = parameters.system_basis_parameters
parameter = "R_int"
para_start = 1
para_end = 300
increment_step = 20

# run_multiple_parameter.run_multiple_parameter(system_basis_parameters, parameter, para_start, para_end, increment_step)
all_spectra = plotting.read_multiple_output(system_basis_parameters, parameter, para_start, para_end, increment_step)

 # %%
def plot_difference_ratio(all_spectra, n_start, n_end, parameters, para, para_start, para_end, increment_step): 
    '''
    Plot energy levels from n_start to n_end from each spectrum 
    difference ratio between numerical and analytical energy over analytical value (of hydrogen atom spectrum) is plotted

    Args:
        all_spectra (list of lists): list of lists containing all the energy spectra
        n_start (int): starting energy level to be plotted (n_start >= 0, inclusive)
        n_end (int): end of range of energy levels to be plotted (inclusive)
    '''

    # Create 2 figures and axis
    fig2, ax2 = plt.subplots(layout="constrained")


    # create a list containing range of a parameter
    parameter_axis = []
    for curr_para in range(para_start, para_end+increment_step, increment_step):
        parameter_axis.append(curr_para)

    # plot energy-para for a range of energy levels
    for i in range(n_start, n_end+1 ,1):
        # difference
        difference_ratios= []
        for spectrum in all_spectra:
            difference_ratios.append((spectrum[i]+0.5/(i+1)**2)/(0.5/(i+1)**2))

        # plot difference_ratio - parameter range
        ax2.plot(parameter_axis, difference_ratios, 'o', label=f"{i}th energy differerence_ratio")

     # Set labels and title for ax2
    ax2.set_xlabel(f"{para}")
    ax2.set_ylabel('difference ratio')
    ax2.set_title('energy difference_ratio - parameter')

    # Display a legend
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    # Add parameters notation on the side
    # convert your dictionary to a string
    params_str = "\n".join(f"{k}: {v}" for k, v in parameters.items())
    fig2.text(2.2, 1, params_str, fontsize=12, ha='left', va='top', transform=ax2.transAxes)
    
    # add a grid background 
    ax2.grid(True)

    # set y-axis to log scale 
    plt.yscale('log')

    # save figure to directory
    file_name2 = f"{n_start}-{n_end}-{ECS12.parametersToString(parameters)}-{para}-{para_start}-{para_end}-{increment_step}-difference_ratio.png"
    # fig2.tight_layout()
    fig2.savefig(f'/users/stud/linzhonglin/linzhonglin/plots/{file_name2}', bbox_inches='tight')
    print(f'figure 2 saved to: /users/stud/linzhonglin/linzhonglin/plots/{file_name2}')

    # Display the plot
    # plt.show(block=False)

n_start_value = 5 # inclusive
n_end_value = 10 # inclusive

plotting.plot_numerical_analytical(all_spectra, n_start_value, n_end_value, system_basis_parameters, parameter, para_start, para_end, increment_step)
plot_difference_ratio(all_spectra, n_start_value, n_end_value, system_basis_parameters, parameter, para_start, para_end, increment_step)

#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline
import argparse
import os.path
import os
import seaborn as sns
import pandas as pd
import fortranformat as ff
from UnitConversion import *
from LaserExcitation import *
from File_reader import *
from Templates import *
import time
import json
from parameters import system_basis_parameters

#########################################################################################
#########################################################################################
#########################################################################################
def parametersToString(parameters):
        '''
        (helper function for createInputFile) 
        Convert parameters dictionary to a string for use in file names
        This string uniquely specifies a input file to the FORTRAN procedure

        Args: 
            parameters (dictionary): see parameter_varied.py
        '''
        temp_dict = parameters.copy()
        temp_dict.pop("N_dig")
        json_string = json.dumps(temp_dict)

        # get rid of special characters for linux filename restrictions
        file_name = json_string.replace(':', '-')
        file_name = file_name.replace("{","")
        file_name = file_name.replace("}","")
        file_name = file_name.replace(",","_")
        file_name = file_name.replace('"','')
        file_name = file_name.replace(" ","")

        return file_name

def createInputFile(parameters):
    '''
    Create the tempate text file to be fed into FORTRAN code.
    Input file is saved to directory specified here in source code. 
    
    Args:: 
        parameters (dict): a dictionary of all the relevant parameters
        >>> exampple:
            system_basis_parameters = {
            # basis
            "NameOfBase": 'Bspline',
            "x_max": 600,
            "N_x": 450,
            "N_y": 30,
            "k_x": 10,
            "k_y": 10,
            # calculation specification
            "m_min": 0,
            "m_max": 0,
            # system 
            "R_int": 1,
            "N_dig": 3
        }
    Returns: 
        None
    '''
    # save input file to the following directory
    directory = '/users/stud/linzhonglin/AMO_TOOLS_USER/ELECTR_STR/B_DAM_ECS/B_DAM_ECS_1e/input/H_H/' 
    
    # give the input file a unique identifier
    # append prefix for H2 system, and suffix for FORTRAN procedure requirement
    filename = ''.join(['H_H_',parametersToString(parameters),'.bde1'])
    fullPath = ''.join((directory,filename))

    
    # calls function bde() defined in Templates.py, which creates a input file for the H2 ion code with the given parameters
    content_exp = bde1(float(parameters['R_int']),parameters['m_min'],parameters['m_max'],parameters['x_max'],parameters['N_x'],parameters['N_y'],parameters['k_x'],parameters['k_y'])
    f_exp1 = open(fullPath,'w')
    f_exp1.write(content_exp)
    f_exp1.close()

    print('created:',fullPath)
    print("finished task c")


#########################################################################################
######################### run calculation ##################################################
#########################################################################################
def check_run_history(parameters):
        '''
        check whether the output file for this input parameters exist in the output direcory 
        return True if exist 
        return False if not exist

        Args: 
            parameters (dictionary): see parameters.py -> sys_basis_parameters
        '''
        output_directory = "/users/stud/linzhonglin/AMO_TOOLS_USER/ELECTR_STR/B_DAM_ECS/B_DAM_ECS_1e/out/H_H/"
        output_file_name = "H_H_" + parametersToString(parameters)
        output_full_path = os.path.join(output_directory, output_file_name)
        if os.path.exists(output_full_path):
            # new found issue: sometimes the file is created but no results inside, so must check # of files inside
            num_files =  len([f for f in os.listdir(output_full_path) if os.path.isfile(os.path.join(output_full_path, f))])
            if num_files == parameters["m_max"]-parameters["m_min"]+2: 
                print(f"'{output_full_path}' \n above otuput file exists, input file ran previously, skip calculation")
                return True
        else:
            print(f"'{output_full_path}' \n input file NOT ran previously")
            return False
        
def runCalculation(parameters):
    '''
    Feed the input file into FORTRAN program for calculation.

    Caution: 
        - this checks for whether same "system_basis_parameters" to FORRTAN has been run before and skip if so
        - this assumes that whatever input file parameters not in "system_basis_parameters" have been kept the same

    Args: 
        parameters (dict): a dictionary of all the relevant parameters
        >>> exampple:
            system_basis_parameters = {
            # basis
            "NameOfBase": 'Bspline',
            "x_max": 300,
            "N_x": 150,
            "N_y": 30,
            "k_x": 10,
            "k_y": 10,
            # calculation specification
            "m_min": 0,
            "m_max": 0,
            # system 
            "R_int": 1,
            "N_dig": 3
        }
    '''
    # first check whether the same input file has been run before 
    # if not check_run_history(parameters): # temporarily commented out to run program with 2nd proton switched on
    program_start_time = time.time()

    # directory where FORTRAN procedure is at
    pathcode='/users/stud/linzhonglin/AMO_TOOLS_USER/ELECTR_STR/B_DAM_ECS/'
    # FORTRAN procedure will append prefix and suffix to find the input file 
    input_file = parametersToString(parameters)

    # run the FORTRAN procedure (assuming already compiled) through a shell command
    ecs1e_task=' '.join((pathcode+'B_DAM_ECS_1e/RUN_B_DAM_ECS_1e.csh','H','H',input_file))
    print(ecs1e_task)
    os.system(ecs1e_task) 
    print("finished task r")

    program_end_time = time.time()
    runtime = program_end_time - program_start_time
    print(f"Runtime: {runtime} seconds")


#########################################################################################
######################### extract info from output ##################################################
#########################################################################################
def readEnergy(parameters, m_show=0, parity_show=""):
    '''
    Read the unformatted output file from FORTRAN procedure 
    and returns EnT -> a list of calculated energy (degeneracy not removed yet)

    Caution: only work for m = 0 now

    Args: 
        parameters (dict): a dictionary of all the relevant parameters
        >>> exampple: 
            system_basis_parameters = {
            # basis
            "NameOfBase": 'Bspline',
            "x_max": 600,
            "N_x": 450,
            "N_y": 30,
            "k_x": 10,
            "k_y": 10,
            # calculation specification
            "m_min": 0,
            "m_max": 0,
            # system 
            "R_int": 1,
            "N_dig": 3
        }
    Optional Args: 
        m_show (int): m quantum number of the state to be read
        parity_show (str): symmetry of the state to be read
            'g': symmetic system -> (“Symmetry of the potential” in input set to 1)
            'u': asymmetric system -> (“Symmetry of the potential” in input set to 2)
    '''

    # output file name is specified in FORTRAN procedure
    # folder that contains output file depend on input file name
    partial_directory='/users/stud/linzhonglin/AMO_TOOLS_USER/ELECTR_STR/B_DAM_ECS/B_DAM_ECS_1e/out/H_H/H_H_'
    directory = partial_directory + parametersToString(parameters)
    
    # depending on whether the system has inversion symmetry
    if parity_show: 
        file_name = "m__" + str(m_show) + "_" + parity_show + ".out1e"
    else: 
        file_name = "m__" + str(m_show) + ".out1e"

    
    full_path = directory + '/' + file_name

    # read the fortran unformatted output file
    print("showing file:"+full_path)
    fout = FortranFile(full_path,'r')

    n_calc_state=fout.read_ints(np.int32)[0]
    print('n_calc_state:', n_calc_state)

    EnT=[0]*n_calc_state
    CxT=[0]*n_calc_state
    CyT=[0]*n_calc_state
    for n in range(n_calc_state):
        EnT[n]=fout.read_reals(float)[0] #,why is this a complex value? second value small
        CxT[n]=fout.read_reals(float)    # Length B_N_X -1  (last 0ne is zero) 
        CyT[n]=fout.read_reals(float)     # Length B_N_y -0   
        
        #since csr stps can be run separatly,
        #this if block is to ensure that when basis parameters are changed, 
        #r step should be run again
        if n==0:
            N_x_f,N_y_f=len(CxT[n])+1,len(CyT[n])
            print('N_x_f,N_y_f(from fortran output):', N_x_f,N_y_f)
            if (N_x_f != parameters['N_x']) or (N_y_f != parameters['N_y']):
                print('ERROR!!!! MISMATCH IN BASE')
                sys.exit()
    fout.close()
    
    return EnT
    
#########################################################################################
######################### plot numerical against numerical #################################
##################################################################################  
def plotEnergy_numerical_analytical(parameters, m_show=0, parity_show="", n_start_value = 5, n_end_value = 10):
    '''
    plots numerical against analytical enery level values on the same plot to quickly check convergence

    Args: 
        parameters (dict): a dictionary of all the relevant parameters
        >>> exampple:
            system_basis_parameters = {
            # basis
            "NameOfBase": 'Bspline',
            "x_max": 600,
            "N_x": 300,
            "N_y": 30,
            "k_x": 10,
            "k_y": 10,
            # calculation specification
            "m_min": 0,
            "m_max": 0,
            # system 
            "R_int": 1,
            "N_dig": 3
        }
        
    optional Args: 
        m_show (int): m quantum number of the state to be read
        parity_show (str): symmetry of the state to be read
            'g': symmetic system -> (“Symmetry of the potential” in input set to 1)
            'u': asymmetric system -> (“Symmetry of the potential” in input set to 2)
        n_start_value (int): starting n value to be plotted (inclusive)
        n_end_value (int): ending n value to be plotted (inclusive)
    '''
    # read output
    EnT = readEnergy(parameters, m_show, parity_show)

    # creat data list to be plotted 
    analytical_energies = []
    numerical_energies = []
    for i in range(n_start_value, n_end_value+1 ,1):
        analytical_energies.append(get_energy_levels(EnT))
        numerical_energies.append(-0.5/(i+1)**2)
    dummy_x_axi = [1] * len(numerical_energies)

    # configure the plot, generate figure
    fig, ax = plt.subplots(layout="constrained")
    for value in analytical_energies: 
        ax.axhline(y=value)
    ax.plot(dummy_x_axi, analytical_energies, linestyle='-', label=f"{i}th analytical energy")

    plt.xlim(0.98, 1.02)
    # plt.ylim(-0.2, 1)

    ax.set_xlabel("dummy axis")
    ax.set_ylabel('energy (atomic unit)')
    ax.set_title('numerical & analytical energy')
    
    plt.grid()
    plt.show()

#########################################################################################
######################### plot wavefunction ##################################################
#########################################################################################
def plotWaveFunctions(parameters, m_show=0, parity_show="", nn=1, hm_range=5, hm_points=70 ):
    '''
    - Read the FORTRAN output of the parameters (unformatted)
    - plot the wave function probability (not normalised) as heat map in Cartesian coordinate.
    - by defaults plots the ground state
    
    TODO: not completed
    
    Args: 
        parameters (dict): a dictionary of all the relevant parameters
        >>> exampple:
            system_basis_parameters = {
            # basis
            "NameOfBase": 'Bspline',
            "x_max": 600,
            "N_x": 300,
            "N_y": 30,
            "k_x": 10,
            "k_y": 10,
            # calculation specification
            "m_min": 0,
            "m_max": 0,
            # system 
            "R_int": 1,
            "N_dig": 3
        }
        
    optional Args: 
        m_show (int): m quantum number of the state to be read
        parity_show (str): symmetry of the state to be read
            'g': symmetic system -> (“Symmetry of the potential” in input set to 1)
            'u': asymmetric system -> (“Symmetry of the potential” in input set to 2)
        nn (int): index of state to be plotted, range? 
        hm_range (int): heat map plot range
        hm_points (int): heat map plot resolution
    '''
    # output file name is specified in FORTRAN procedure
    # folder that contains output file depend on input file name
    partial_directory='/users/stud/linzhonglin/AMO_TOOLS_USER/ELECTR_STR/B_DAM_ECS/B_DAM_ECS_1e/out/H_H/H_H_'
    directory = partial_directory + parametersToString(parameters)
    
    # depending on whether the system has inversion symmetry
    if parity_show: 
        file_name = "m__" + str(m_show) + "_" + parity_show + ".out1e"
    else: 
        file_name = "m__" + str(m_show) + ".out1e"
    
    full_path = directory + '/' + file_name

    # read the fortran unformatted output file
    print("showing file:"+full_path)
    fout = FortranFile(full_path,'r')

    n_calc_state=fout.read_ints(np.int32)[0]
    print('n_calc_state:', n_calc_state)

    EnT=[0]*n_calc_state
    CxT=[0]*n_calc_state
    CyT=[0]*n_calc_state
    for n in range(n_calc_state):
        EnT[n]=fout.read_reals(float)[0] #,why is this a complex value? second value small,type: list of int
        CxT[n]=fout.read_reals(float)    # Length B_N_X -1  (last 0ne is zero), type: list of lists 
        CyT[n]=fout.read_reals(float)     # Length B_N_y -0, type: list of lists   
        
        #since csr stps can be run separatly,
        #this if block is to ensure that when basis parameters are changed, 
        #r step should be run again
        if n==0:
            N_x_f,N_y_f=len(CxT[n])+1,len(CyT[n])
            print('N_x_f,N_y_f(from fortran output):', N_x_f,N_y_f)
            if (N_x_f != parameters['N_x']) or (N_y_f != parameters['N_y']):
                print('ERROR!!!! MISMATCH IN BASE')
                sys.exit()
    fout.close()



    def symmetry_dependent_x(x, lambda_, parameters, C_coeffs):
        '''
        symmetry dependent xi function (2.35) in Yulian's thesis

        Args: 
            x: xi coordinate in spheroidal coordinate system
            lambda_ (int): angular momentum quantum number, n = 0,1,2,3... -> needed to determine symmetry
            parameters (dict): the parameters dictionary as per usual -> see parameters.py
            C_coeffs (list): list of coefficients from FOTRAN output 
        '''

        # knots
        knots_x = np.concatenate((
            [1]*(parameters['k_x']-1),
            np.linspace(1,parameters['x_max'],
                        parameters['N_x']-parameters['k_x']+2), 
                        [parameters['x_max']]*(parameters['k_x']-1)
                                  ))
        
        # degree of BSpline 
        # -> because scipy.interpolate.BSpline() takes degree as input rather than order
        kx = parameters['k_x']-1
        
        # adapt FORTRAN output coefficients list for input to scipy.interpolate.BSpline()
        Coeffs_adapted = np.concatenate((C_coeffs,[0])) # because FORTRAN output enforces last coeff to be 0

        print(f"Calling BSpline() for xi coordinate axis, with parameters: \n\
        number of knots: {len(knots_x)} \n\
        number of spline coefficients: {len(C_coeffs)} \n\
            B-spline degree: {kx}")
        
        output = (x**2-1)**(lambda_/2) * BSpline(knots_x, Coeffs_adapted, kx)(x)
        
        return output
    
    # test plotting Figure 2.2 in Yulian's thesis
    def plot_symmetry_dependent_x(lambda_, parameters, C_coeffs):
        '''
        plots the symmetry_dependent_x for a range of values specified in source code

        Args: 
            lambda_ (int): related to absolute value of m quantum number of state plotted
            parameters (dict): see parameters.py
            C_coeffs (np.ndarray): FORTRAN output coefficients
        
        '''
        # create data to be plotted
        x_values = np.linspace(1, 50, 10000)
        y_values = [symmetry_dependent_x(x, lambda_, parameters, C_coeffs) for x in x_values]
        fig, ax = plt.subplots()
        ax.plot(x_values, y_values)
        
        plt.title('Symmetry Dependent Function xi coordinate')
        plt.xlabel('xi coodinate')
        plt.ylabel('symmetry_dependent_x')
        plt.title('Plot of symmetry_dependent_x')
        plt.grid(True)
        # set x-axis to log scale 
        plt.xscale('log')
        # plt.xlim(0, 50)
        # plt.ylim(-0.2, 1)
        plt.show()


    def symmetry_dependent_y(y, lambda_, parity, parameters, C_coeffs):
        '''
        symmetry dependent xi function (2.31) in Yulian's thesis

        Args: 
            x: xi coordinate in spheroidal coordinate system
            lambda_ (int): angular momentum quantum number, n = 0,1,2,3... -> needed to determine symmetry
            parity (str): specify state symmetry
                'g': weierstrass p = 1
                'u': weierstrass p = -1
            parameters (dict): the parameters dictionary as per usual -> see parameters.py
            C_coeffs (list): list of coefficients from FOTRAN output 
        '''
        parameters = system_basis_parameters

        # knots 
        knots_y=np.concatenate((
            [-1]*(parameters['k_y']-1),
            np.linspace(-1,1,2*parameters['N_y']-parameters['k_y']+2), 
            [1]*(parameters['k_y']-1)
            ))
        
        # degree of BSpline 
        # -> because scipy.interpolate.BSpline() takes degree as input rather than order
        ky = parameters['k_y']-1

    
        # adapt FORTRAN output coefficients list for input to scipy.interpolate.BSpline()
        Coeffs_adapted_first_half = np.concatenate((C_coeffs, np.zeros_like(C_coeffs)))
        Coeffs_adapted_second_half =  np.concatenate((np.zeros_like(C_coeffs), C_coeffs[::-1]))

        print(f"Calling BSpline() for eta coordinate axis, with parameters: \n\
              number of knots: {len(knots_y)} \n\
                len(C_coeffs): {len(C_coeffs)} \n\
                len(Coeffs_adapted_first_half): {len(Coeffs_adapted_first_half)} \n\
                len(Coeffs_adapted_second_half): {len(Coeffs_adapted_second_half)} \n\
                    B-spline degree: {ky}")    
        
        if parity == 'g': 
            p = 1
        elif parity == 'u': 
            p = -1
        output = (1-y**2)**(lambda_/2) * (BSpline(knots_y, Coeffs_adapted_first_half, ky)(y) + (-1)**lambda_ * p * BSpline(knots_y, Coeffs_adapted_second_half, ky)(y))
        
        return output

    # test plotting Figure 2.2 in Yulian's thesis
    def plot_symmetry_dependent_y(lambda_, parity, parameters, C_coeffs):
        '''
        plots the symmetry_dependent_y funtion for a range of values specified in source code 

        Args: 
            lambda_ (int): related to absolute value of m quantum number of state plotted
            parameters (dict): see parameters.py
            C_coeffs (np.ndarray): FORTRAN output coefficients
        
        '''
        # create data to be plotted
        x_values = np.linspace(-1, 1, 10000)
        y_values = [symmetry_dependent_y(x, lambda_, parity, parameters, C_coeffs) for x in x_values]
        fig, ax = plt.subplots()
        ax.plot(x_values, y_values)
        
        plt.title('Symmetry Dependent Function eta coordinate')
        plt.xlabel('eta coordinate')
        plt.ylabel('symmetry_dependent_y')
        plt.title('Plot of symmetry_dependent_y')
        plt.grid(True)
        # set x-axis to log scale 
        # plt.xscale('log')
        # plt.xlim(0, 50)
        # plt.ylim(-0.2, 1)
        plt.show()

    lambda_ = 0
    n = 1
    C_coeffs_x = CxT[n]
    C_coeffs_y = CyT[n]
    plot_symmetry_dependent_x(lambda_, system_basis_parameters, C_coeffs_x)
    # plot_symmetry_dependent_y(lambda_, 'u', system_basis_parameters, C_coeffs_y)
    
    def wavefunction(nn, lambda_, parameters, C_coeffs, x1, y1):
        '''
        Returns the wave function value at Cartesian coordinate position x1, y1
        for the number nn state
        (NOT normalised)

        Args: 
            nn (int): index of state in FORTRAN output to be calcualted
            lambda_ (int): angular momentum quantum number, n = 0,1,2,3... -> needed to determine symmetry   
            C_coeffs (list): list of coefficients from FOTRAN output          
            x1 (int): x coordinate axis in Cartesian coordinate system, range?
            y1 (int): y coordinate axis in Cartesian coordinate system, range? 
        '''
        xi,yi=cart2sphr(x1,y1,parameters['R_int'])
        phi = symmetry_dependent_x(xi, lambda_, parameters, C_coeffs) * symmetry_dependent_y(yi, )
        return phi

    def plot_probability_hm(parameters, hm_range, hm_points):
        '''
        plot the (not normalised) probability distribution as a heat map

        Args: 
            parameters (dict): to get R_int
            hm_range (int): ploting range -> center to side distance of a square
            hm_points (int): determines plot resolution
        '''

        fig, axs= plt.subplots(1,parameters['m_max']+1)

        # create heat map grid
        hm_x1 = np.linspace(-hm_range*parameters['R_int'], hm_range*parameters['R_int'], hm_points)
        hm_y1 = np.linspace(-hm_range*parameters['R_int'], hm_range*parameters['R_int'], hm_points)
        M = np.zeros((len(hm_x1),len(hm_y1)))

        # add values to heat map grid
        for y1 in range(len(hm_y1)):
            for x1 in range(len(hm_x1)):
                print(f"Now calculating value at \n\
                    Cartesian coordinate: x={x1}, y={y1} \n\
                    State index: {nn} \n")
                M[y1,x1]= wavefunction(nn, hm_x1[x1],hm_y1[y1])**2
        
        # plot heat map 
        # sns.heatmap(M, cmap='Blues',xticklabels=False,yticklabels=False,cbar=False,ax=axs[0])
        plt.show()

    
    print("finished task w")

'''
    def bsp_x(x_coord,c):    #input:b-spline parameters of length B_N-2
            c0=np.concatenate((c,[0])) # because FORTRAN output enforces last coeff to be 0
            bsp = BSpline(knots_x, c0, parameters['k_x']-1)
            return np.array(bsp(x_coord))  #output: y-values  
    
    def bsp_y(y_coord,c,gu):    #input:b-spline parameters of length B_N-2
        
        if gu=='g':
            c0=np.concatenate((c,c[::-1]))  # was guessed 
        if gu=='u':
            c0=np.concatenate((-c,c[::-1]))  # was guessed 


        bsp = BSpline(knots_y, c0, parameters['k_y']-1)
        return np.array(bsp(y_coord))  #output: y-values 

    def elec(z_,r_):
                xi,yi=cart2sphr(z_,r_,parameters['R_int'])
                phi=bsp_x(xi,CxT[parameters['nn']])*bsp_y(yi,CyT[parameters['nn']],parameters['parity_show'])
                return phi
    

    fig, axs= plt.subplots(2,parameters['m_max']+1)

    hm_range=5 #plot range
    hm_points=70 #plot resolution

    hm_zz = np.linspace(-hm_range*parameters['R_int'], hm_range*parameters['R_int'], hm_points)
    hm_rr = np.linspace(-hm_range*parameters['R_int'], hm_range*parameters['R_int'], hm_points)

    M=np.zeros((len(hm_zz),len(hm_rr)))
    for r in range(len(hm_rr)):
        for z in range(len(hm_zz)):
            M[r,z]= elec(hm_zz[z],hm_rr[r])**2

    if parameters['parity_show']=='g':        
        sns.heatmap(M, cmap='Blues',xticklabels=False,yticklabels=False,cbar=False,ax=axs[0])
    if parameters['parity_show']=='u':       
        sns.heatmap(M, cmap='Greens',xticklabels=False,yticklabels=False,cbar=False,ax=axs[0])

    plt.show()
    print("finished task w")
'''


#########################################################################################
######################### helper functions ##################################################
#########################################################################################
def get_energy_levels(EnT):
            '''
            returns a list of the energy levels by considering 1st value of each degenerate energy level

            Args: 
                EnT (list): list of calculated energies
            '''
            def nth_energy(EnT, n):
                '''
                Returns the first value as nth energy for each degenerate level
                If index out of range, return False

                Caution: for now, this code only caters to m=0 case -> i.e. gorund state energy starts from 0.5
                
                TODO:
                    - make this work for all m values
                    - upper bound for n not coded here

                Args: 
                    EnT (list): list of calcualted energy from FORTRAN otuput
                    n (int): [1, upper bound?]
                '''
                index = int((n*(n-1))/2)
                if n <= 0: 
                    raise Exception("n value error: n >= 1")
                if index < len(EnT): 
                    return EnT[index]
                else: 
                    return False
                
            # consider 1st value of each degenerate level 
            spectrum = []
            for n in range(1,len(EnT)): 
                if nth_energy(EnT, n):
                    spectrum.append(nth_energy(EnT,n))
            return spectrum

plotEnergy_numerical_analytical(parameters, m_show=0)
# %%
