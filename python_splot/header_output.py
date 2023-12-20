"""
Takes the stored data from composition_jpt_entropy_reader and writes
it all to an output file called sphToMesa.out which can be used
for further analysis of the stars.

I recommend using Jacky Tran's Hypermongo tool which is a versatile
graphing tool that can be used to plot all variables from SPH
against each other very nicely.

Jacky's email: jtran9148@gmail.com
Jacky's GitHub: https://github.com/IndigenousAlien/Hypermongo
Jacky's Senior Comp Paper: Modeling Tidal Disruption Encounters and Analysis with a Data Visualization Tool - Hypermongo
     The paper can be found on https://dspace.allegheny.edu/
     (It may need to be requested in order to obtain access)


Charles Gibson
Allegheny College & Northwestern University
Department of Physics & CIERA
Illinois Space Grant Consortium
08/26/2023
"""

import numpy as np
from composition_jpt_entropy_reader import entropy_reader

def header_output(interp_data,star_num=1):
    
    xq=interp_data['xq']
    q=interp_data['q']
    r=interp_data['r']
    rho=interp_data['rho']
    P=interp_data['P']
    A=interp_data['A']
    T=interp_data['T']
    jrot=interp_data['jrot']
    H1=interp_data['H1']
    He3=interp_data['He3']
    He4=interp_data['He4']
    C12=interp_data['C12']
    N14=interp_data['N14']
    O16=interp_data['O16']
    Ne20=interp_data['Ne20']
    Mg24=interp_data['Mg24']

    print(f'WRITING TO sphToMesa{star_num}.out')
    
    with open(f'sphToMesa{star_num}.out','w') as f:
        f.write('xq     q     r     rho     P     A     T     jrot     H1     He3     He4     C12     N14     O16     Ne20     Mg24\n')

        for i in range(len(xq)):
            f.write(f'{xq[i]}  {q[i]}  {r[i]}  {rho[i]}  {P[i]}  {A[i]}  {T[i]}  {jrot[i]}  {H1[i]}  {He3[i]}  {He4[i]}  {C12[i]}  {N14[i]}  {O16[i]}  {Ne20[i]}  {Mg24[i]}\n')

