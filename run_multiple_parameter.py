import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import time
import argparse
import ECS12
import parameters

def run_multiple_parameter(parameters, para, para_start, para_end, increment_step):
    '''
    Run multiple calculations for a particulat parameter with the range of parameter values specified

    Args: 
        - para (string): parameter to be varied across a range
    '''
    start_time = time.time()
    for curr_para in range(para_start, para_end+increment_step, increment_step):
        # break loop if exceeding range
        if curr_para > para_end+increment_step*2: 
            break

        print(f"Now running calculation for {para} = {curr_para}")

        # modify parameter
        parameters[para] = curr_para
        
        ## create input file and run calulation 
        ECS12.createInputFile(parameters)
        ECS12.runCalculation(parameters)

        print(f"Finished running calculation for {para} = {curr_para}")

    print("All calcualtions done")
    end_time = time.time()
    print(f"Total runtime: {end_time - start_time} seconds")

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description=
                                     'This python script runs the FORTRAN procedure \
                                     for a range of value for a specified parameter \n \
                                     Caution: edit other basis parameters in the source code')

    # Add the arguments
    parser.add_argument('parameter',
                        choices=['R_int', 'N_x', 'x_max'],
                        type=str,
                        help='The parameter to be varied. Parameters that can be varied: %(choices)s.')

    parser.add_argument('para_start',
                        type=int,
                        help='The start of the varied range. (inclusive)')

    parser.add_argument('para_end',
                        type=int,
                        help='The end of the varied range. (inclusive)')
    
    parser.add_argument('increment_step',
                        type=int,
                        help='Step by which the parameter is incremented')

    # Execute the parse_args() method
    args = parser.parse_args()


    system_basis_parameters = parameters.system_basis_parameters

    run_multiple_parameter(system_basis_parameters, args.parameter, args.para_start, args.para_end, args.increment_step)

if __name__ == '__main__': 
    main()