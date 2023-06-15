"""
This code uses the data read in from the readit.py file as well as the
splines from composition_spline.py to determine the composition of each particle.
It then writes all values to composition.sph.

This will be used as reference for future interactions with the star, as we can
assume that the composition of each particle will remain the same between interactions.

Charles Gibson
Allegheny College
Department of Physics
05/26/2023
"""

import numpy as np

def composition_fit(data, comp_splines):

    # necessary imported data for next calculations
    ntot = data['ntot']
    rho = data['rho']
    am = data['am']
    h1 = comp_splines['h1']
    he3 = comp_splines['he3']
    he4 = comp_splines['he4']
    c12 = comp_splines['c12']
    n14 = comp_splines['n14']
    o16 = comp_splines['o16']
    ne20 = comp_splines['ne20']
    mg24 = comp_splines['mg24']

    nbins = 1601 # number of bins for sorting data
    nbinsbf = 100 # number of bins for the best fit data: higher values means higher resolution
    alogrhomin = -7.9 # 10^-10 is the lowest density
    alogrhomax = 2.1 # 10^10 is the highest density
    amasstot = np.sum(am)

    # initializes amrho
    amrho = np.zeros(nbins)
    
    for i in range(ntot):
        alogrhoi = np.log10(rho[i]) # finds the log base 10 of each density value
        # determines which bin to put the particle
        nbini = int((alogrhoi-alogrhomin) * (nbins-1.0) / (alogrhomax - alogrhomin) + 1.0)

        # checks if the bin is outside of the range of possible bins
        if(nbini > nbins or nbini <  1):
            print("PROBLEM 3")

        # these are leftover codes from the original Fortran. This method also works but is
        # much less efficient
        # this line adds the mass element for the particle to the bin that the particle falls in
        # for k in range(1, nbini):
        #     amrho[k] += am[i]

        amrho[nbini] += 0.5 * am[i] # adds 1/2 the mass of the particle to its corresponding bin

    mycompfile = open("composition.sph", "w")

    # adds the particle mass to all the bins before the corresponding bin to account for mass enclosed
    # in isodense surfaces
    for i in range(nbins):
        for k in range(0, i):
            amrho[k] += amrho[i]*2

    for i in range(ntot):
        alogrhoi = np.log10(rho[i]) # finds the log of the density of each particle
        nbini = (alogrhoi - alogrhomin) * (nbins - 1) / (alogrhomax - alogrhomin) + 1 # index of the bin that the particle is in
        amonmrho = amrho[int(nbini)] / amasstot # mass fraction of the particle based on the bin it falls into
        mycompfile.write(f'{h1(amonmrho)}    {he3(amonmrho)}    {he4(amonmrho)}    {c12(amonmrho)}\
    {n14(amonmrho)}    {o16(amonmrho)}    {ne20(amonmrho)}    {mg24(amonmrho)}\n')
