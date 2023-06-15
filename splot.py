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
user = str(getpass.getuser()) # gets the username of the user (example: vanderzyden01) to use in line 4 (for the splot folder location)
sys.path.insert(0, f'/home/{user}/splot_directories/python_splot') # this location may change for the future, if so this line must be edited accordingly
# these files must all be imported after the above two lines
from bestfit import bestfit
from grapefree_readit import grapefree_readit
from readit import readit
from readit_collision import readit_collision
from compsph_readit import compsph_readit
from composition_jpt_entropy_reader import entropy_reader
from read_orig_profile_comp import composition_reader
from composition_spline import comp_spline
from composition_fit import composition_fit
from component_reader import component_reader
from bound_particles import bound_particle_data
from component_best3 import compbest3

def main():
    ### Similar to pplot.f
    # print to the user to ask for the option number for bestfit, composition, etc.
    option = int(input("Which option would you like to do?\n    Generate dat files for MESA [1] \n\
    Generate composition.sph [2] \n    Analyze a star from a collision [3] \n\
    Determine components [4] \n: "))
    if option == 1 or option == 3:
        mode = input("Which type of input file would you like to feed to MESA?\n    DT\n    PT\n    DE\n    DP\n: ")
        neos = int(input("Which type of eos is your simulation using?\n    Analytic (default SPH) [1]\n    MESA [2]\n: "))
    nnit_input_start = int(input("Starting out file (int): "))
    nnit_input_end = int(input("Ending out file (int, must be >= the start value): "))
    frequency = int(input("Step size between files: "))
    nnit_input = nnit_input_start
    while nnit_input <= nnit_input_end:
        if nnit_input_end >= nnit_input_start:
            print(f"--------------------------Output file {nnit_input}--------------------------")
            readit_data = readit_collision(nnit_input, 4) # iform is always 4 (usually inputted from a separate routine, but I have it set to 4)
            if option == 1:
                comp_data = compsph_readit()
                bestfit(readit_data, comp_data, neos)
                entropy_reader(mode)
            elif option == 2:
                q, elements = composition_reader()
                splines = comp_spline(q, elements)
                composition_fit(readit_data, splines)
            elif option == 3:
                # Option 2 still needs to be run with the original relaxation to create composition.sph

                component_data = component_reader(nnit_input) # makes a list of the values of the component for each particle
                orig_comp_data = compsph_readit() # reads the composition data from composition.sph - move outside the loop
                # finds the bound particle and composition data to be passed to bestfit_total.py
                bound_data, bound_composition_data = bound_particle_data(readit_data, component_data, orig_comp_data)
                # creates composition.dat, entropy.dat, and angular_momentum.dat
                bestfit(bound_data, bound_composition_data, neos)
                entropy_reader(mode)
            elif option == 4:
                compbest3(nnit_input,readit_data)
            print("Output data analysis completed.")
            nnit_input += frequency
            if nnit_input_end < nnit_input_start:
                print("End value was smaller than start value. Please try again.")
    ###

if __name__ == "__main__":
    main()
