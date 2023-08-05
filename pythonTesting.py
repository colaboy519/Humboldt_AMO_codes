#%%
import ECS12
from parameters import system_basis_parameters
parameters = system_basis_parameters
from scipy.io import FortranFile
from UnitConversion import p_str
import numpy as np
import sys

m_show = 0

# first check parameters have been run 
if ECS12.check_run_history(parameters):
    # output file name is specified in FORTRAN procedure
    # folder that contains output file depend on input file name
    partial_directory='/users/stud/linzhonglin/AMO_TOOLS_USER/ELECTR_STR/B_DAM_ECS/B_DAM_ECS_1e/out/H_H/H_H_'
    directory = partial_directory + ECS12.parametersToString(parameters)

    # depending on symmetry of the problem: comment out one of the following lines accordingly 
    # symmetic system
    # file_name = "m__"+str(parameters['m_show'])+"_"+parameters['parity_show']+".out1e"
    # asymmetric system
    file_name = "m__"+str(m_show)+".out1e"

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

else:
    print("parameters not run yet")
# %%
import numpy as np

n = 5
i = 3

arrary = (np.arange(n) == i).astype(float)
print(arrary)
# %%
import numpy as np

def double_size_with_zeros(arr):
    zeros = np.zeros_like(arr)
    new_arr = np.concatenate((arr, zeros))
    return new_arr

# Example usage
my_array = np.array([1, 2, 3, 4])
doubled_array = double_size_with_zeros(my_array)
print("double array: ", doubled_array)
print("my_array: ", my_array)


# %%
