"""
This code uses the data read in from the out****.sph file as well as
the splines / linear interpolation from composition_interpolator.py
to determine the composition of each particle.
It then writes all values to sph.composition.

This will be used as reference for future interactions with the star,
as we can assume that the composition of each particle will remain the
same between interactions.

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

    nbins = 16001 # number of bins for sorting data
    # number of bins for the best fit data
    nbinsbf = 100 # higher values means higher resolution
    alogrhomin = -13 # 10^-13 is the lowest density
    alogrhomax = 5 # 10^5 is the highest density
    amasstot = np.sum(am)

    # initializes amrho
    amrho = np.zeros(nbins)
    
    for i in range(ntot):
        # finds the log base 10 of each density value
        alogrhoi = np.log10(rho[i])
        # determines which bin to put the particle
        nbini = int(
            (alogrhoi-alogrhomin) * (nbins-1.0)
            / (alogrhomax - alogrhomin) + 1.0)

        # checks if the bin is outside of the range of possible bins
        if(nbini > nbins or nbini <  1):
            print("PROBLEM 3")

        # adds 1/2 the mass of the particle to its corresponding bin
        amrho[nbini] += 0.5 * am[i]

    mycompfile = open("sph.composition", "w")

    # adds the particle mass to all the bins before the corresponding
    # bin to account for mass enclosed in isodense surfaces
    for i in range(nbins):
        for k in range(0, i):
            amrho[k] += amrho[i]*2

    for i in range(ntot):
        # finds the log of the density of each particle
        alogrhoi = np.log10(rho[i])
        # index of the bin that the particle is in
        nbini = (alogrhoi - alogrhomin) * (nbins - 1) \
            / (alogrhomax - alogrhomin) + 1
        # mass fraction of the particle based on the bin it falls into
        amonmrho = amrho[int(nbini)] / amasstot
        mycompfile.write(
f'{h1(amonmrho)}    {he3(amonmrho)}    {he4(amonmrho)}    {c12(amonmrho)}  \
  {n14(amonmrho)}    {o16(amonmrho)}    {ne20(amonmrho)}    {mg24(amonmrho)}\n')
