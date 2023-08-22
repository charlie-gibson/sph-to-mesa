"""Hans Vanderzyden Python3 Splot File for Analysis of StarSmasher (Main)
The purpose of this file is to be a main file that calls in other subroutines (in /data/user/splot folder created by user) needed to analyze the output data from StarSmasher.
Currently only reads in the data (additional functionality is in the process)

---------------------------------------------------------------------------------------------------
Updated to read in data and do data analysis to write MESA friendly files
for future MESA evolution from an SPH simulation. There is an option to run this
post processing routine (1) and an option to determine the composition of each
particle (2) which can be used on a relaxed star for future evaluation of other interactions

Charles Gibson
Allegheny College
Department of Physics
"""

import sys
import getpass
import os
import glob
import re
import numpy as np
user = str(getpass.getuser()) # gets the username of the user (example: vanderzyden01) to use in line 4 (for the splot folder location)
sys.path.insert(0, f'/home/{user}/sph-to-mesa/python_splot/') # this location may change for the future, if so this line must be edited accordingly
# these files must all be imported after the above two lines
from bestfit_smoothed2 import bestfit
from readit_collision import readit_collision
from compsph_readit import compsph_readit
from composition_jpt_entropy_reader import entropy_reader
from read_orig_profile_comp import composition_reader
from composition_spline import comp_spline
from composition_fit import composition_fit
from component_reader import component_reader
from bound_particles import bound_particle_data
from component_best3 import compbest3
from pa_plot import pa_plot
from energy_graph import v

def main():
    ### Similar to pplot.f
    # print to the user to ask for the option number for bestfit, composition, etc.
    option = int(input("Which option would you like to do?\n    Generate MESA files - no collision [1] \n\
    Generate composition.sph [2] \n    Generate MESA files - collision [3] \n\
    Determine components [4] \n    pa plot [5] \n    v plot [6] \n: "))
    if option == 1 or option == 3:
        mode = input("Which type of input file would you like to feed to MESA?\n    DT\n    PT\n    DE\n    DP\n: ")
        neos = int(input("Which type of eos is your simulation using?\n    Analytic (default SPH) [1]\n    MESA [2]\n: "))
        comp_data = compsph_readit()
    if option == 2 or option ==5:
        profile_num=int(input('What profile is used in the SPH relaxation (just number)? '))
    if option==3:
        component_val=int(input("Which star are you analyzing? "))
    if option == 6:
        v()
        raise SystemExit
    nnit_input_start = int(input("Starting out file (int): "))
    nnit_input_end = int(input("Ending out file (int, must be >= the start value): "))
    frequency = int(input("Step size between files: "))
    nnit_input = nnit_input_start
    while nnit_input <= nnit_input_end:
        if nnit_input_end >= nnit_input_start:
            print(f"--------------------------Output file {nnit_input}--------------------------")
            readit_data = readit_collision(nnit_input, 4) # iform is always 4 (usually inputted from a separate routine, but I have it set to 4)
            if option == 1:
                bestfit(readit_data, comp_data, neos)
                entropy_reader(mode)
            elif option == 2:
                q, elements = composition_reader(profile_num)
                splines = comp_spline(q, elements)
                composition_fit(readit_data, splines)
            elif option == 3:
                # Option 2 still needs to be run with the original relaxation to create composition.sph
                component_data = component_reader(nnit_input) # makes a list of the values of the component for each particle
                # finds the bound particle and composition data to be passed to bestfit_total.py
                bound_data, bound_composition_data = bound_particle_data(readit_data, component_data, comp_data,component_val)
                # creates composition.dat, entropy.dat, and angular_momentum.dat
                bestfit(bound_data, bound_composition_data, neos)
                entropy_reader(mode)
            elif option == 4:
                try:
                    file_list=glob.glob('comp****.sph')
                    #print(file_list)
                    largest_number=0
                    largest_file='comp'
                    for file in file_list:
                        #print(file)
                        file_number=int(file[4:8])
                        #print(file_number)
                        if file_number >= largest_number and file_number <= nnit_input:
                            largest_number=int(np.ceil((file_number+0.1)/10)*10)
                            largest_file=file
                    #print(largest_file)
                    with open(largest_file,'r') as f:
                        icomp=[]
                        for line in f:
                            vals=line.split()
                            icomp.append(vals[3])
                    #print(icomp)
                except FileNotFoundError:
                    print('NO EXISTING comp****.sph FILE')
                    n1=readit_data['n1']
                    icomp1=[1]*n1
                    n2=readit_data['n2']
                    icomp2=[2]*n2
                    icomp=icomp1+icomp2
                    largest_number=0
                except KeyError:
                    print('NO EXISITNG comp****.sph FILE\nONE STAR IN SPH SIMULATION')
                    ntot=readit_data['ntot']
                    icomp=[1]*ntot
                    largest_number=0
                except Exception as e:
                    print('An error occurred: ',e)
                print('icomp calculated')
                print('STARTING OUTPUT FILE: ',largest_number)
                dtout=readit_data['dtout']
                if dtout>10:
                    step=int(dtout)
                else:
                    step=int(10/dtout)
                for i in range(largest_number,nnit_input+1,step):
                    write=False
                    print(f'\n--------------------------------------Output File {i}----------------------------------------------')
                    readit_data=readit_collision(i,4)
                    if i==nnit_input or i%50==0:
                        write=True
                    icomp=compbest3(i,readit_data,icomp,write)
                    if i+step - nnit_input > step:
                        write=True
                        icomp=compbest3(nnit_input,readit_data,write)
            elif option == 5:
                pa_plot(readit_data,profile_num,nnit_input)
            print("Output data analysis completed.")
            nnit_input += frequency
            if nnit_input_end < nnit_input_start:
                print("End value was smaller than start value. Please try again.")
    ###

if __name__ == "__main__":
    main()
