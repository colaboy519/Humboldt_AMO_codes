# Humboldt_AMO_codes
This repository stores the python codes used to run the FORTRAN calculation codes for hydrogen, hydrogen molecule in AMO group at Humboldt. 

Since these codes are more like helper codes to make use of the FORRTAN codes, they don't contain any sensitive calculation-related information. Hence I am putting them up here in public for ease of reference and whoever finds them useful in case someone is following up with my project done Humbolt, AMO group in summer 2023 related to simulting Rydberg atom orbitals in an optical lattice using numerical methods.

# Baisc instructions on using the python scripts
1. Specify the basis and system parameters in parameters.py
2. For calculation of a single set of parameters, run ECS12.py
3. For calculation of a range of a particular parameter, run run_multiple_parameter.py 
    1. To plot eigen values against range of parameters for all calculation outputs, run plotting.py


## Code devleopment protocol
1. To save calculation time, each calculation output is given a unique name to prevent redundant calculations

## User-interface scripts
### ECS12.py
This python script can separately create input file for FORTRAN procedure, feed input file into running FORTRAN procedure, and read output from FORTRAN procedure.

### Hatom.py
Caution: this script is under development!
equivalent of ECS12.py for hydrogen atom system code rather than hydogen molecule system code. As the two FORTRAN procedures takes in difference input template and outputs with difference formates, the input and output reading of ECS12.py need to be adapted for this to work.

### run_multiple_parameter.py
This python script runs the FORTRAN procedure for a range of value for a specified parameter.
Caution: edit other basis parameters in the source code

### plotting.py
-> This python script is used after run_multiple_parameter.py 
This script plots the energy spectra of the range of parameters on the same plot plots are saved with file name with format: {n_start}-{n_end}-{ECS12.parametersToString(parameters)}-{para}-{para_start}-{para_end}-{increment_step}-difference.png\

### H2ionConversionSearch.py
Caution: this script is under development!
This script was meant to automatically search for minimum basis for the H2 ion system that is able to satisfy a certain specified convergence criterion. 

## Helper Python Scripts
### parameters.py
Helper python script that contains the basis and system parameters that can be changed before the python script feeds the input file into FORTRAN procedure. 
-> example: 
    system_basis_parameters = {
            # basis
            "NameOfBase": 'Bspline',
            "x_max": 1000,
            "N_x": 500,
            "N_y": 30,
            "k_x": 10, # order of xi coordinate BSplines
            "k_y": 10,
            # calculation specification
            "m_min": 0,
            "m_max": 0,
            # system 
            "R_int": 141,
            "N_dig": 3
        }

### Templates.py
Helper python script that contains the input template with keys in parameters.py as variable.

### UnitConversion.py
Helper python script for unit conversion and coordinate transformation between Cartesian and Elliptical coordinates. 
