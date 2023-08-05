import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from File_reader import *
import ECS12
import argparse
import parameters

def read_multiple_output(parameters, para, para_start, para_end, increment_step):
    '''
    Read and store energy levels from the range of a parameter
    returns a list of list

    Caution: only work ofor m=0 now (single output file in output folder)

    Args: 
        - para (string): parameter to be varied across a range
    '''
    all_spectra = []
    for curr_para in range(para_start, para_end+increment_step, increment_step):
        print(f"Now reading output for {para} = {curr_para}")

        # modify parameter
        parameters[para] = curr_para

        # add energy spectrum from this set of parameters of output file can be found
        if ECS12.check_run_history(parameters):
            EnT = ECS12.readEnergy(parameters,0)
            energies = ECS12.get_energy_levels(EnT)
            all_spectra.append(energies)
        else: 
            output_directory = "/users/stud/linzhonglin/AMO_TOOLS_USER/ELECTR_STR/B_DAM_ECS/B_DAM_ECS_1e/out/H_H/"
            output_file_name = "H_H_" + ECS12.parametersToString(parameters)
            output_full_path = os.path.join(output_directory, output_file_name)
            raise Exception(f"{output_full_path} \n the above input file has not been calculated")
    print("Finished reading all energy spectra.")
    return all_spectra
    

def plot_numerical_analytical(all_spectra, n_start, n_end, parameters, para, para_start, para_end, increment_step): 
    '''
    Plot energy levels from n_start to n_end from each spectrum 
    Numerical is plotted together with analytical value (of hydrogen atom spectrum)

    Args:
        all_spectra (list of lists): list of lists containing all the energy spectra
        n_start (int): starting energy level to be plotted (n_start >= 0, inclusive)
        n_end (int): end of range of energy levels to be plotted (inclusive)
    '''

    # Create 2 figures and axis
    fig1, ax1 = plt.subplots(layout="constrained")

    # create a list containing range of a parameter
    parameter_axis = []
    for curr_para in range(para_start, para_end+increment_step, increment_step):
        parameter_axis.append(curr_para)

    # plot energy-para for a range of energy levels
    for i in range(n_start, n_end+1 ,1):
        # break loop before exceeding range (because mind not have been calcualted)
        if i > n_end+1: 
            break

        # ith numerical energy 
        numerical_energies = []
        # ith analytical energy
        analytical_energies = []
        for spectrum in all_spectra:
            numerical_energies.append(spectrum[i])
            analytical_energies.append(-0.5/(i+1)**2)
        
        # plot numerical & analytical energy - parameter range
        ax1.plot(parameter_axis, numerical_energies,'o', label=f"{i}th numerical energy ")
        ax1.plot(parameter_axis, analytical_energies, linestyle='-', label=f"{i}th analytical energy")

    # Set labels and title for ax1
    ax1.set_xlabel(f"{para}")
    ax1.set_ylabel('energy (atomic unit)')
    ax1.set_title('numerical & analytical energy - parameter')

    # Display a legend
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    # Add parameters notation on the side
    # convert your dictionary to a string
    params_str = "\n".join(f"{k}: {v}" for k, v in parameters.items())
    fig1.text(1.8, 1, params_str, fontsize=12, ha='left', va='top', transform=ax1.transAxes)
    
    # add a grid background 
    ax1.grid(True)

    # save figures to directory
    file_name1 = f"{n_start}-{n_end}-{ECS12.parametersToString(parameters)}-{para}-{para_start}-{para_end}-{increment_step}.png"
    # fig1.tight_layout()
    fig1.savefig(f'/users/stud/linzhonglin/linzhonglin/plots/{file_name1}', bbox_inches='tight')
    print(f'figure 1 saved to: /users/stud/linzhonglin/linzhonglin/plots/{file_name1}')

    # Display the plot
    # plt.show(block=False)

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
    fig2.text(1.8, 1, params_str, fontsize=12, ha='left', va='top', transform=ax2.transAxes)
    
    # add a grid background 
    ax2.grid(True)

    # set y-axis to log scale 
    # plt.yscale('log')

    # save figure to directory
    file_name2 = f"{n_start}-{n_end}-{ECS12.parametersToString(parameters)}-{para}-{para_start}-{para_end}-{increment_step}-difference_ratio.png"
    # fig2.tight_layout()
    fig2.savefig(f'/users/stud/linzhonglin/linzhonglin/plots/{file_name2}', bbox_inches='tight')
    print(f'figure 2 saved to: /users/stud/linzhonglin/linzhonglin/plots/{file_name2}')

    # Display the plot
    # plt.show(block=False)

def main(): 
    # Create the parser
    parser = argparse.ArgumentParser(description=\
                                    'This python script is used after run_multiple_parameter.py \n\
                                        This script plots the energy spectra of the range of parameters on the same plot \n\
                                            plots are saved with file name with format: {n_start}-{n_end}-{ECS12.parametersToString(parameters)}-{para}-{para_start}-{para_end}-{increment_step}-difference.png')
    
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

    parser.add_argument('-n_s','--n_start_value', default=0,
                        type=int,
                        help='starting energy level to be plotted, n_start>=0')
    
    parser.add_argument('-n_e','--n_end_value', default=5,
                        type=int,
                        help='end of energy level range to be plotted (n_end >= 0, inclusive) \n\
                              n_end exceeding number of energy levels calcualted will raise error')

    # Execute the parse_args() method
    args = parser.parse_args()


    system_basis_parameters = parameters.system_basis_parameters

    all_spectra = read_multiple_output(system_basis_parameters, args.parameter, args.para_start, args.para_end, args.increment_step)
    plot_numerical_analytical(all_spectra, args.n_start_value, args.n_end_value, system_basis_parameters, args.parameter, args.para_start, args.para_end, args.increment_step)
    plot_difference_ratio(all_spectra, args.n_start_value, args.n_end_value, system_basis_parameters, args.parameter, args.para_start, args.para_end, args.increment_step)

if __name__ == '__main__': 
    main()