#!/usr/bin/python3

"""
Hans Vanderzyden
Python3 Splot File for Analysis of StarSmasher (Main)
The purpose of this file is to be a main file that calls in other
subroutines(in /data/user/splot folder created by user)
needed to analyze the output data from StarSmasher.
Currently only reads in the data
(additional functionality is in the process)

-----------------------------------------------------------------------
This is the main hub for StarSmasher analysis

Updated to allow additional functionality from other data analysis
codes. Options for analysis include:
0  - Writes an output file in ASCII format
1  - Reconstructs a stellar model, ignoring smoothing of composition
     and angular momentum
2  - Writes an sph.composition file based on a relaxation
3  - Reconstructs a stellar model, incorporating smoothing of
     composition and angular momentum
4  - Evaluates the bound state of each particle
5  - Rewrites MESA files with a given thermodynamic pair, especially
     for the entropy?.dat file, using a pre-existing bestfit?.sph file
6  - Generates a particle plot for relaxations, similar to pa in smx
7  - Generates an energy plot by reading in the energy?.sph file
8  - Forces a model to be in hydrostatic equilibrium
9  - Calculates the orbital eccentricity and relative v_infty
10 - Uses existing *.dat files and writes all MESA input files
     including inlists for a MESA construction
11 - Combines options 2, 4, 3, and 10 for a relaxation
12 - Combines options 4, 3, and 10 for a collision
13 - Writes a file with particle data consistent with how
     StarSmasher writes data

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
user = str(getpass.getuser()) # gets username of user
# path where sph-to-mesa exists
# sys.path.insert(0, f'/jet/home/{user}/sph-to-mesa/python_splot/')
sys.path.insert(0, f'home/{user}/sph-to-mesa/python_splot')

# these and all new files must be imported after the above two lines
from ascii_output import ascii_output
from bestfit import bestfit
from readit_collision import readit_collision
from compsph_readit import compsph_readit
from composition_jpt_entropy_reader import entropy_reader
from read_orig_profile_comp import composition_reader
# from composition_spline import comp_spline
from composition_interpolator import comp_spline
from composition_fit import composition_fit
from component_reader import component_reader
from bound_particles import bound_particle_data
from component_best3 import compbest3
from pa_plot import pa_plot
from energy_graph import v
from input_reader import sph_input_reader
from neighbors import neighbors
from header_output import header_output
from hydrostatic_equilibrium_neos1 import hse_func
from eccentricity import eccentricity
from write_mesa_inlist import write_mesa_inlist
from write_bound_particles import write_bound_particles
from write_bound_composition import write_bound_composition

# path where sph-to-mesa lives
path='/home/{user}/sph-to-mesa/python_splot'

def option0(nnit_input, readit_data, header_data):

    """
    writes SPH data to ASCII format
    """

    ascii_output(nnit_input, readit_data, header_data)

def option1(readit_data, sph_input):

    """
    Reconstructs stellar model without smoothing any values
    such as composition and angular momentum
    """

    # entropy pairs for entropy?.dat
    mode = input("Which type of input file would you like to feed to MESA?\n\
    Entropy    		                 [S]  \n\
    Density-Temperature			 [DT] \n\
    Pressure-Temperature      	         [PT] \n\
    Density-Specific Internal Energy     [DE] \n\
    Density-Pressure (Experimental)	 [DP] \n\
: ")

    # evaluates particle composition from sph.composition
    comp_data = compsph_readit()

    # determines which equation of state to use
    try:
        neos = sph_input['neos']
    except:
        neos = 1

    # writes bestfit?.sph and resolved, interpolated data
    bestfit(readit_data, comp_data, neos, sph_input)
    header_data=entropy_reader(mode,path)
    header_output(header_data)

def option2(readit_data, sph_input):

    """
    Writes an sph.composition file
    """

    profile_num = sph_input['profilefile']
    q, elements = composition_reader(profile_num)
    splines = comp_spline(q, elements, spline=False)
    composition_fit(readit_data, splines)

def option3(nnit_input, readit_data, sph_input, mode=-1, component_val=-1):

    """
    Reconstructs the stellar model while smoothing composition
    and angular momentum according to the Wendland C4 kernel
    """

    if mode == -1:
        mode = input(\
"Which type of input file would you like to feed to MESA?\n\
    Entropy    		                 [S]  \n\
    Density-Temperature			 [DT] \n\
    Pressure-Temperature      	         [PT] \n\
    Density-Specific Internal Energy     [DE] \n\
    Density-Pressure (Experimental)	 [DP] \n\
: ")

    if component_val == -1:
        component_val=int(input("Which star are you analyzing? "))

    comp_data = compsph_readit()

    try:
        neos = sph_input['neos']
    except:
        neos = 1

    component_data = component_reader(nnit_input)
    composition_data, jrot_data, temp_data, p_data = neighbors(
        readit_data,comp_data,sph_input,component_val,component_data,neos
    )
    print('FINDING BOUND PARTICLES')
    bound_data, bound_composition_data, jrot_data, temp_data, p_data = \
        bound_particle_data(
            readit_data, component_data, composition_data, jrot_data,
            temp_data, p_data, component_val
        )
    print('ENTERING BESTFIT')
    bestfit(
        bound_data, bound_composition_data, neos, sph_input, jrot_data,
        temp_data, p_data, component_val)
    header_data = entropy_reader(mode, path, component_val)
    header_output(header_data, component_val)

    return bound_data, bound_composition_data

def option4(
        nnit_input, nnit_input_end, readit_data, icomp,
        write_interval, internal_energy_fraction
    ):

    """
    Determines the bound components of each particle employing the
    approach presented in Lombardi, et al. (2006)

    Select 0 for internal energy contribution for only mechanical
    (gravitational potential + kinetic) energy contributions as in
    Nandez, et al. (2014)
    """

    icomp_old = icomp
    write = False
    print(f'  READING OUTPUT FILE {nnit_input}')
    if nnit_input % write_interval == 0 or nnit_input == nnit_input_end:
        write = True
    
    icomp, convergence = compbest3(
        nnit_input, readit_data, icomp_old, write, internal_energy_fraction
    )
    if not convergence:
        icomp = icomp_old
        n += 1
        if n >= 5:
            print('CANNOT CONVERGE\nEXITING CODE')
            raise SystemExit

    return icomp

def option5():

    """
    Uses an existing bestfit?.sph file to construct a new entropy?.dat
    file with a new thermodynamic pair

    Entropy, S, only works when using an eos file with entropy, setting
    neos=2 in sph.input
    """

    mode = input("Which type of input file would you like to feed to MESA?\n\
    Entropy    		                 [S]  \n\
    Density-Temperature			 [DT] \n\
    Pressure-Temperature      	         [PT] \n\
    Density-Specific Internal Energy     [DE] \n\
    Density-Pressure (Experimental)	 [DP] \n\
: ")

    component_val=int(input("Which star are you analyzing? "))

    header_data = entropy_reader(mode, path, component_val)
    header_output(header_data, component_val)

def option6(readit_data, sph_input, nnit_input):

    """
    Creates a plot comparing particle data to the profile used for
    stellar relaxation

    Similar to using supermongo and creating pa plots with smx
    """

    profile_num=sph_input['profilefile']
    pa_plot(readit_data, profile_num, nnit_input)

def option7():

    """
    Reads energy?.sph and plots everything in it
    including all energies, angular momentum, and entropy
    """

    num = int(input("Energy File Number: "))
    v(num)
    raise SystemExit

def option8(mode=-1, component_val=-1):

    """
    Forces a star into hydrostatic equilibrium by solving
    the equation dP/dr = -g*rho

    Accounts for slight rotational effects well below breakup
    rotation
    """

    if mode == -1:
        mode = input("Which type of input file would you like to feed to MESA?\n\
    Entropy    		                 [S]  \n\
    Density-Temperature			 [DT] \n\
    Pressure-Temperature      	         [PT] \n\
    Density-Specific Internal Energy     [DE] \n\
    Density-Pressure (Experimental)	 [DP] \n\
: ")
    if component_val == -1:
        component_val=int(input("Which star are you analyzing? "))

    hse_func(component_val)
    header_data=entropy_reader(mode, path, component_val, 1)
    header_output(header_data, component_val, 1)

def option9(readit_data, nnit_input):

    """
    Uses orbital mechanics to return eccentricity and
    relative velocity at infinity

    Equations come from Chapter 8 of
    Classical Dynamics of Particles and Systems
    Thornton & Marion (5th Edition)
    """

    eccentricity(readit_data, nnit_input)

def option10(
        nnit_input, readit_data, sph_input, bound_data=None,
        bound_composition_data=None, mode=-1, component_val=-1, hse_input=-1
    ):

    """
    Writes simple MESA inlists that allow for stellar evolution
    of relaxations in MESA

    Inlists are written for mesa-r24.03.1 and may need to be
    updated for future releases
    """

    if mode == -1:
        mode = input(
"Which type of input file would you like to feed to MESA?\n\
    Entropy    		                 [S]  \n\
    Density-Temperature			 [DT] \n\
    Pressure-Temperature      	         [PT] \n\
    Density-Specific Internal Energy     [DE] \n\
    Density-Pressure (Experimental)	 [DP] \n\
: ")

    if component_val == -1:
        component_val=int(input("Which star are you analyzing? "))

    if hse_input == -1:
        hse_input = int(input('HSE?\n yes [1]\n no  [0]\n: '))
    if hse_input == 0:
        hse = False
    elif hse_input == 1:
        hse = True

    try:
        neos = sph_input['neos']
    except:
        neos = 1
        
    comp_data = compsph_readit()

    component_data = component_reader(nnit_input)
    if not bound_data:
        composition_data, jrot_data, temp_data, p_data = neighbors(
            readit_data, comp_data, sph_input, component_val,
            component_data, neos)
        print('FINDING BOUND PARTICLES')
        bound_data, bound_composition_data, jrot_data, temp_data, p_data = \
            bound_particle_data(
                readit_data, component_data, composition_data, jrot_data,
                temp_data, p_data, component_val
            )
    write_mesa_inlist(
        bound_data, bound_composition_data, mode, component_val, hse
    )
    raise SystemExit

def option11(nnit_input, nnit_input_end, readit_data, sph_input):

    """
    Runs a complete stellar reconstruction of the final relaxed
    star in StarSmasher

    Creates an sph.composition file,
    generates entropy1.dat, angular_momentum1.dat,
    and composition1.dat using smoothed data,
    and writes MESA inlist files
    """

    mode = input("Which type of input file would you like to feed to MESA?\n\
    Entropy    		                 [S]  \n\
    Density-Temperature			 [DT] \n\
    Pressure-Temperature      	         [PT] \n\
    Density-Specific Internal Energy     [DE] \n\
    Density-Pressure (Experimental)	 [DP] \n\
: ")

    component_val=int(input("Which star are you analyzing? "))

    icomp = readit_data['ntot'] * [1]

    option2(readit_data, sph_input)

    option4(
        nnit_input, nnit_input_end,readit_data, icomp,
        write_interval=nnit_input_end, internal_energy_fraction=0
    )
    bound_data, bound_composition_data = option3(
        nnit_input, readit_data, sph_input, mode, component_val
    )
    option10(
        nnit_input, readit_data, sph_input, bound_data,
        bound_composition_data, mode, component_val, hse_input=0
    )

    raise SystemExit

def option12(
        nnit_input, readit_data, sph_input, icomp,
        mode=[-1,-1], component_val=[1, 2]
    ):

    """
    Does a complete stellar reconstruction of a collision product
    in StarSmasher

    Creates an sph.composition file from
    an sph.composition1 and sph.composition2 file (if sph.composition
    does not yet exist),
    evaluates bound components of particles,
    generates entropy1.dat, angular_momentum1.dat,
    and composition1.dat using smoothed data,
    and writes MESA inlist files
    """

    for comp_val in component_val:

        # print('STAR',component_num)
        print('MODE',mode[comp_val-1])

        component_num = 0

        # determines how many particles are bound to star?
        for i in icomp:
            if i == comp_val:
                component_num += 1

        if component_num == 0:
            print(f'NO PARTICLES BOUND TO STAR {comp_val}')
            continue

        if mode[comp_val - 1] == -1:

            mode = input(
"Which type of input file would you like to feed to MESA?\n\
    Entropy    		                 [S]  \n\
    Density-Temperature			 [DT] \n\
    Pressure-Temperature      	         [PT] \n\
    Density-Specific Internal Energy     [DE] \n\
    Density-Pressure (Experimental)	 [DP] \n\
: ")

        bound_data, bound_composition_data = option3(
            nnit_input, readit_data, sph_input, mode[comp_val-1], comp_val
        )
        option10(
            nnit_input, readit_data, sph_input, bound_data,
            bound_composition_data, mode[comp_val-1], comp_val, hse_input=0
        )
        os.system(f'mv MESA_input_files MESA_input_files{comp_val}')

def option13(nnit_input, header_data, readit_data, component_val=-1):

    """
    Writes a new bound*.sph file with the same format as an out*.sph
    which can be used as a restartrad file or an sph.start?u for
    another StarSmasher run
    """

    if component_val == -1:
        component_val=int(input("Which star are you analyzing? "))
    component_data = component_reader(nnit_input)
    composition_data = compsph_readit()
    bound_data, bound_composition = bound_particle_data(
        data=readit_data, component=component_data,
        composition=composition_data, component_val=component_val,
        all_data=False
    )
    write_bound_particles(nnit_input, header_data, bound_data, component_val)
    write_bound_composition(nnit_input, bound_composition, component_val)    

def main():

    ### based on pplot.f

    path=f'/home/{user}/python_splot'
    option = int(input("Which option would you like to do?\n\
    Write ascii output                         [0] \n\
    Generate MESA files - no smoothing         [1] \n\
    Generate sph.composition                   [2] \n\
    Generate MESA files - smoothing            [3] \n\
    Determine components                       [4] \n\
    Generate MESA files - existing bestfit     [5] \n\
    pa plot                                    [6] \n\
    v plot                                     [7] \n\
    Force Hydrostatic Equilibrium              [8] \n\
    Calculate Eccentricity                     [9] \n\
    Write MESA inlist                          [10] \n\
    Complete relaxation reconstruction         [11] \n\
    Complete collision reconstruction          [12] \n\
    Write outfile with only bound particles    [13] \n: "))

    icomp = [-1]
    write_interval = -1
    internal_energy_fraction = -1
    complete = False

    if option != 12:
        nnit_input_start = int(input("Starting out file: "))
        nnit_input_end = int(input("Ending out file: "))
        frequency = int(input("Step size: "))

    elif option == 12:
        # gets last output file
        cwd = os.getcwd()
        files = os.listdir(cwd)
        outfiles = [
            file for file in files if
            file.startswith('out') and file.endswith('.sph')
        ]
        outnums = []
        for out in outfiles:
            outnums.append(int(out.replace('out','').replace('.sph','')))
        # automatcially decides to use the first and last output file
        # in the collision for the calculation
        nnit_input_start = min(outnums)
        nnit_input_end = max(outnums)
        frequency = 1
        mode = []
        for i in range(2):
            modei = input(f"MESA ENTROPY FILE FOR STAR {i+1}\n\
    Entropy    		                 [S]  \n\
    Density-Temperature			 [DT] \n\
    Pressure-Temperature      	         [PT] \n\
    Density-Specific Internal Energy     [DE] \n\
    Density-Pressure (Experimental)	 [DP] \n\
: ")
            mode.append(modei)

    nnit_input = nnit_input_start

    # used to ensure correct naming of the file
    if nnit_input < 10:
        fname = f'out000{nnit_input}.sph'
    elif nnit_input < 100:
        fname = f'out00{nnit_input}.sph'
    elif nnit_input < 1000:
        fname = f'out0{nnit_input}.sph'
    else:
        fname = f'out{nnit_input}.sph'

    print('------------------- sph.input ------------------')
    sph_input=sph_input_reader()
    print('\n\n')
    try:
        neos=sph_input['neos']
    except:
        neos=1
    
    orig_option = option

    if option == 12:
        # writes sph.composition using sph.composition1
        # and sph.composition2 is sph.composition does not exist
        files = os.listdir()
        if 'sph.composition' not in files:
            os.system(
                'cat sph.composition1 sph.composition2 >> sph.composition'
            )
        # automatically runs option 4 to create all comp****.sph files
        option = 4
        write_interval = 1000 # decrease to get more comp****.sph files
        internal_energy_fraction = 0

    while nnit_input <= nnit_input_end:
        if nnit_input_end >= nnit_input_start:
            print(
f"--------------------------Output file {nnit_input}-------------------------"
            )
           
            if option !=5:
                if option == 0:
                    readit_data, header_data = readit_collision(
                        fname, 4, True
                    )
                elif option == 13:
                    readit_data, header_data = readit_collision(
                        fname, 3, True
                    )
                else:
                    readit_data = readit_collision(fname, 4)

            if option == 0:
                option0(nnit_input, readit_data, header_data)

            if option == 1:
                option1(readit_data, sph_input)
            
            if option == 2:
                option2(readit_data, sph_input)

            if option == 3:
                option3(nnit_input, readit_data, sph_input)

            if option == 4:

                if write_interval == -1:
                    write_interval = int(input(
"Frequency of writing comp*.sph files: "
                    ))
                if internal_energy_fraction == -1:
                    internal_energy_fraction = float(input(
"Enter fraction of internal energy to use when determining bound mass: "
                    ))
            
                if icomp[0] == -1:
                    # begins by listing all comp****.sph files
                    # in the directory
                    try:
                        file_list=glob.glob('comp****.sph')
                        # largest_number=0
                        largest_number = nnit_input_start
                        if file in file_list:
                            for file in file_list:
                                print(file_list)
                                # gets number of each output file
                                file_number = int(str(
                                    file.replace('comp','')
                                    ).replace('.sph',''))
                                print(file_number, largest_number, nnit_input)
                                if (file_number >= largest_number and
                                    file_number < nnit_input):
                                    largest_number = file_number
                                    largest_file = file
                                elif (file_number >= largest_number and
                                      file_number > nnit_input):
                                    largest_number = file_number
                                    if largest_number > 9999:
                                        largest_file = \
                                            f'comp{largest_number:05d}.sph'
                                    else:
                                        largest_file = \
                                            f'comp{largest_number:04d}.sph'
                                elif (file_number == largest_number and
                                      file_number == nnit_input):
                                    largest_number = file_number
                                    largest_file = file
                    except:
                        largest_number = nnit_input_start
                        largest_file = f'comp{nnit_input_start:04d}.sph'

                    print(largest_file)

                    # allows the component determination to begin from
                    # the most recent comp****.sph file
                    if nnit_input != largest_number:
                            # used to ensure correct naming of the file
                        if largest_number < 10:
                            fname = f'out000{largest_number}.sph'
                        elif largest_number < 100:
                            fname = f'out00{largest_number}.sph'
                        elif largest_number < 1000:
                            fname = f'out0{largest_number}.sph'
                        else:
                            fname = f'out{largest_number}.sph'

                        readit_data = readit_collision(fname, 4)
                    nnit_input = largest_number

                    # tries to open comp****.sph file 
                    try:
                        with open(largest_file,'r') as f:
                            icomp=[]
                            n=0
                            for line in f:
                                if n==0:
                                    n+=1
                                else:
                                    vals=line.split()
                                    icomp.append(int(vals[3]))
                        print(np.where(icomp==1))
                        print(np.where(icomp==2))
                    except:
                        # if there are no comp****.sph files for a
                        # 2+ body interaction
                        try:
                            print('NO EXISTING comp****.sph FILE')
                            n1=readit_data['n1']
                            icomp1=[1]*n1
                            n2=readit_data['n2']
                            icomp2=[2]*n2
                            icomp=icomp1+icomp2
                            largest_number=0
                        except KeyError:
                            # if there are no comp****.sph files
                            # for a stellar relaxation
                            print(
'NO EXISITNG comp****.sph FILE\nONE STAR IN SPH SIMULATION'
                            )
                            ntot=readit_data['ntot']
                            icomp=[1]*ntot
                            largest_number=0

                        print('ICOMP CALCULATED')
                    print('STARTING OUTPUT FILE: ',largest_number)

                else:
                    icomp = icomp
                    
                icomp = option4(
                    nnit_input, nnit_input_end, readit_data, icomp,
                    write_interval, internal_energy_fraction
                )

                if orig_option == 12:
                    if nnit_input == nnit_input_end:
                        option = orig_option


            if option == 5:
                option5()

            if option == 6:
                option6(readit_data, sph_input, nnit_input)

            if option == 7:
                option7()

            elif option == 8:
                option8()

            elif option == 9:
                option9(readit_data, nnit_input)

            elif option == 10:
                option10(nnit_input, readit_data, sph_input)
            
            if option == 11:
                option11(nnit_input, nnit_input_end, readit_data, sph_input)

            if option == 12:
                option12(nnit_input, readit_data, sph_input, icomp, mode=mode)

            if option == 13:
                option13(
                    nnit_input, header_data, readit_data, component_val=-1
                )


            if nnit_input == nnit_input_end:
                break
            nnit_input += frequency
            if nnit_input > nnit_input_end and not complete:
                if orig_option == 12:
                    option = orig_option
                    complete = False
                else:
                    complete = True
                    nnit_input = nnit_input_end
            if nnit_input_end < nnit_input_start:
                print(
"End value was smaller than start value. Please try again."
                )
            print('ANALYSIS COMPLETE')

if __name__ == "__main__":
    main()
