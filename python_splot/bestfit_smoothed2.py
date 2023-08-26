"""
This code calls the py file from splot_directories.

Using the read in data, this code sorts every particle based on the density
into many bins. It also uses the MESA EOS table or the analytical EOS equations
to determine the temperature and pressure of each particle before it is binned.
Finally, it writes these data to a best fit table.

The original Fortran file (written by James Lombardi) states:
Analyze the results of a collision run and output important
numbers for a summary table

Charles Gibson & James Lombardi
Allegheny College
Department of Physics
"""

# data is read in using readit.py and passed through splot.py
def bestfit(data, comp_data, neos,sph_input):

    # from readit import readit
    from get_temperature import get_temperature
    import numpy as np
    from eos_func import read_eos
    from eos_func import useeostable

    # uses the lists in data (a dictionary) to determine the variables for calculations
    ntot = int(data['ntot']) # number of particles
    x = data['x'] # x position
    y = data['y'] # y position
    z = data['z'] # z position
    am = data['am'] # mass
    hp = data['hp'] # smoothing length
    rho = data['rho'] # density
    vx = data['vx'] # x velocity
    vy = data['vy'] # y velocity
    vz = data['vz'] # z velocity
    vxdot = data['vxdot'] # x acceleration
    vydot = data['vydot'] # y acceleration
    vzdot = data['vzdot'] # z acceleration
    u = data['u'] # specific internal energy
    udot = data['udot'] # d/dt[u]
    grpot = data['grpot'] # gravitational potential
    meanmolecular = data['meanmolecular'] # mean molecular weight
    cc = data['cc']
    divv = data['divv']
    # aa = data['aa']
    # bb = data['bb']
    # dd = data['dd']
    # dt = data['dt'] * 0.018445 * 24 * 60 * 60 # to get dt in seconds

    for i, values in enumerate(comp_data):
        if i == 0:
            h1 = values
        elif i == 1:
            he3 = values
        elif i == 2:
            he4 = values
        elif i == 3:
            c12 = values
        elif i == 4:
            n14 = values
        elif i == 5:
            o16 = values
        elif i == 6:
            ne20 = values
        elif i == 7:
            mg24 = values

    print("DATA PASSED TO FUNCTION")
    
    # variables
    nbins = 16001 # number of bins for sorting data
    nbinsbf = 100 # number of bins for the best fit data: higher values means higher resolution
    alogrhomin = -13 # 10^-10 is the lowest density
    alogrhomax = 4 # 10^10 is the highest density

    # constants used for calculations
    gravconst = 6.67390e-08 # m^3 / (kg s^2)
    munit = 1.9891e33 # kg
    runit = 6.9599e10 # m
    tunit = np.sqrt(runit**3/(gravconst * munit)) # s
    boltz = 1.380658e-16 # erg/kelvin
    crad = 2.997924580e10 # cm/sec  NOTE: card has a different meaning in MESA
    planck = 6.6260755e-27 # gram cm^2/sec
    crad2 = crad**2 # cm^2/s^2
    sigma = np.pi**2*boltz*(boltz*2*np.pi/planck)**3/60/crad2 #cgs
    arad = 4.0*sigma/crad #cgs
    qconst = 1.5*boltz/arad #cgs

    temp = np.zeros(ntot)
    pgas = np.zeros(ntot)
    prad = np.zeros(ntot)
    ptot = np.zeros(ntot)
    # number is ntot - 1 to account for one particle being a black hole
    # I would like to find a way to make this line more general to avoid having to edit
    # the code for each new run
    rhocgs = np.zeros(ntot)
    ucgs = np.zeros(ntot)
    qval = np.zeros(ntot)
    rval = np.zeros(ntot)

    # rhocgs = rho * munit/runit**3
    # ucgs = u * gravconst * munit/runit

    # print(arad)
    # print(rhocgs)
    # print(ucgs)
    # print(qval)
    # print(rval)

    if neos == 1:
        rhocgs = rho*munit/runit**3
        # print(rhocgs)
        for i in range(ntot):
            ucgs[i] = u[i]*gravconst*munit/runit
        # print(ucgs)
        qval = qconst*rhocgs/meanmolecular
        # print(qval)
        rval = -ucgs*rhocgs/arad
        # print(rval)

        for n in range(ntot):
            try:
                temp[n] = get_temperature(qval[n], rval[n])
                pgas[n] = rhocgs[n] * boltz * temp[n] / (meanmolecular[n])
                # prad[n] = arad * temp[n]**4/3
                ptot[n] = pgas[n] # + prad[n]
            except:
                print("ERROR CALCULATING TEMPERATURE OR PRESSURE")

    elif neos == 2:
        eosfile=sph_input['eosfile']
        eos_data = read_eos(eosfile)
        # print(eos_data)
        for i in range(ntot):
            try:
                rhocgs[i] = rho[i]*munit/runit**3
                ucgs[i] = u[i]*gravconst*munit/runit
                # qval[i] = qconst*rhocgs[i]/meanmolecular[i]
                # rval[i] = -ucgs[i]*rhocgs[i]/arad

                temp[i] = useeostable(eos_data=eos_data, ucgs=ucgs[i], rhocgs=rhocgs[i], xxx=h1[i], which=0)
                ptot[i] = useeostable(eos_data=eos_data, ucgs=ucgs[i], rhocgs=rhocgs[i], xxx=h1[i], which=2)
                pgas[i] = ptot[i] - arad * temp[i]**4/3
            except:
                print("problem using eostable", i, rho[i]*munit/runit**3, u[i]*gravconst*munit/runit, useeostable(eos_data=eos_data, ucgs=ucgs[i], rhocgs=rhocgs[i], xxx=h1[i], which=0))
                raise SystemExit
    
    # for i in temp:
    #     print(i)

    # Initialize array where amrho is the total mass in each bin separated by density
    amrho = np.zeros(nbins)

    # initiating the total mass and center of mass positions and velocities
    amasstot = np.sum(am)
    xc = np.sum(x * am) / amasstot
    yc = np.sum(y * am) / amasstot
    zc = np.sum(z * am) / amasstot
    vxc = np.sum(vx * am) / amasstot
    vyc = np.sum(vy * am) / amasstot
    vzc = np.sum(vz * am) / amasstot

    # Bin #1 corresponds to rho = 10**alogrhomin
    # Bin #nbins corresponds to rho = 10**alogrhomax
    # Log(rho) = alogrhomin + (alogrhomax - alogrhomin) * (bin# - 1)/(nbins - 1)

    ammrhoavg = np.zeros(nbinsbf)
    rhoavg = np.zeros(nbinsbf)
    uiavg = np.zeros(nbinsbf)
    ravg = np.zeros(nbinsbf)
    amavg = np.zeros(nbinsbf)

    # loops through every particle
    for i in range(ntot):
        alogrhoi = np.log10(rho[i]) # finds the log base 10 of each density value
        # determines which bin to put the particle
        nbini = int((alogrhoi-alogrhomin) * (nbins-1.0)/(alogrhomax - alogrhomin))

        # checks if the bin is outside of the range of possible bins
        if(nbini >= nbins or nbini <  0):
            print("PROBLEM 3", nbini, alogrhoi)

        # these are leftover codes from the original Fortran. This method also works but is
        # much less efficient
        # this line adds the mass element for the particle to the bin that the particle falls in
        # for k in range(1, nbini):
        #     amrho[k] += am[i]

        amrho[nbini] += 0.5 * am[i] # adds 1/2 the mass of the particle to its corresponding bin

    # adds the particle mass to all the bins before the corresponding bin to account for mass enclosed
    # in isodense surfaces
    for i in range(nbins):
        for k in range(0, i):
            amrho[k] += amrho[i]*2

    print("PARTICLES BINNED")

    # prints the total mass and the center of mass
    print(f"Total Mass = {amasstot}")
    print(f"Center of Mass Position: ({xc}, {yc}, {zc})")
    print(f"Center of Mass Velocity: ({vxc}, {vyc}, {vzc})")
    
    # initializes the lists that will be used for further binning
    ammrhoavg = np.zeros(nbinsbf)
    amavg = np.zeros(nbinsbf)
    uiavg = np.zeros(nbinsbf)
    rhoavg = np.zeros(nbinsbf)
    ravg = np.zeros(nbinsbf)
    jrotavg = np.zeros((nbinsbf, 3))
    radius = 0
    tavg = np.zeros(nbinsbf)
    pavg = np.zeros(nbinsbf)
    r_array = np.zeros(ntot)
    amavgcomp=np.zeros(nbinsbf)

    # composition arrays
    h1avg = np.zeros(nbinsbf)
    he3avg = np.zeros(nbinsbf)
    he4avg = np.zeros(nbinsbf)
    c12avg = np.zeros(nbinsbf)
    n14avg = np.zeros(nbinsbf)
    o16avg = np.zeros(nbinsbf)
    ne20avg = np.zeros(nbinsbf)
    mg24avg = np.zeros(nbinsbf)

    # loops over every particle
    for i in range(ntot):
        alogrhoi = np.log10(rho[i]) # finds the log of the density of each particle
        nbini = int((alogrhoi - alogrhomin) * (nbins - 1) / (alogrhomax - alogrhomin)) # index of the bin that the particle is in
        nbinimin=nbini
        smoothing_factor=5
        for ii in range(nbini-1,-1,-1):
            if amrho[ii] >= min(amrho[nbini] + smoothing_factor * am[i], amrho[0]):
                nbinimin=ii
                break
            nbinimax=nbini
        for ii in range(nbini+1,nbins):
            if amrho[ii] <= max(amrho[nbini] - smoothing_factor * am[i], 0):
                nbinimax=ii
                break
        amonmrho = amrho[nbini] / amasstot # mass fraction of the particle based on the bin it falls into
        index = int(amonmrho * nbinsbf)
        amonmrhomin=amrho[nbinimin]/amasstot # mass fraction of the particle based on the bin it fall sinto
        indexmax=min(int(amonmrhomin*nbinsbf),nbinsbf-1) # Note: "max" on left-hand side because smaller density is larger index value
        amonmrhomax=amrho[nbinimax]/amasstot # mass fraction of the particle based on the bin it falls into
        indexmin=int(amonmrhomax*nbinsbf) # Note:"min" on left-hand side because larger density is smaller index value
        if index >= nbinsbf:
            print("Warning: index >= nbinsbf", index)
            index = nbinsbf - 1
        
        if index < 0:
            print("Warning: index < 0", index)
            raise SystemExit

        # note that all of these multiply by mass but are divided by another mass later on
        amavg[index] += am[i] # adds the mass for each particle to its corresponding index
        uiavg[index] += u[i] * am[i] # this gives the internal energy for each index
        rhoavg[index] += rho[i] * am[i] # this gives the density times mass for each index
        rvec = np.array([x[i] - xc, y[i] - yc, z[i] - zc]) # this provides the r vector of the particle
        vvec = np.array([(vx[i]-vxc), (vy[i]-vyc), (vz[i]-vzc)]) # this gives the velocity vector of the particle
        r = np.linalg.norm(rvec) # this is the length of the r vector from the center of mass
        # print(r)
        r_array[i] += r
        v = np.linalg.norm(vvec) # this is the speed of the v vector relative to the center of mass
        ravg[index] += r * am[i] # this gives the distance times the mass of the particle
        jrotavg[index] += am[i] * np.cross(rvec, vvec) # this gives the angular momentum times the mass of the particle
        tavg[index] += temp[i] * am[i]
        pavg[index] += pgas[i] * am[i]

        # composition data

        midpoint = amonmrho * nbinsbf - 0.5
        max_distance = max(indexmax-midpoint, midpoint-indexmin)+1
        denominator = 0
        #print(indexmin, indexmax)
        for index2 in range(indexmin, indexmax + 1):
            distance = abs(index2 - midpoint)
            q = distance / max_distance
            kernel_value = 2 * q**3 - 3 * q**2 + 1 # A simple function that has a slope of zero at the middle (q=0) and edges
            denominator += kernel_value
        for index2 in range(indexmin, indexmax + 1):
            distance = abs(index2 - midpoint)
            max_distance = indexmax + 1 - midpoint
            q = distance / max_distance
            kernel_value = 2 * q**3 - 3 * q**2 + 1 # A simple function that has a slope of zero at the edges and middle
            fraction = kernel_value / denominator
            #print(len(amavgcomp))
            #print(index2)
            amavgcomp[index2] += am[i]*fraction
            h1avg[index2] += h1[i] * am[i] * fraction
            he3avg[index2] += he3[i] * am[i] * fraction
            he4avg[index2] += he4[i] * am[i] * fraction
            c12avg[index2] += c12[i] * am[i] * fraction
            n14avg[index2] += n14[i] * am[i] * fraction
            o16avg[index2] += o16[i] * am[i] * fraction
            ne20avg[index2] += ne20[i] * am[i] * fraction
            mg24avg[index2] += mg24[i] * am[i] * fraction

        # finds these values for the outermost particle
        if r > radius:
            radius = r
            aimax = u[i]
            rhomin = rho[i]
            jrotmax = am[i] * np.cross(rvec, vvec)
            tmax = temp[i]
            pmax = pgas[i]

            h1max = h1[i]
            he3max = he3[i]
            he4max = he4[i]
            c12max = c12[i]
            n14max = n14[i]
            o16max = o16[i]
            ne20max = ne20[i]
            mg24max = mg24[i]

        ammrhoavg[index] = ammrhoavg[index] + amonmrho * am[i]

    #print(f'H1 average: {h1avg}')
        
    # calculates quantities of values in certain regions of the star

    #print(r_array)
    #print(am)
    #print(radius)
    enclosed_mass = []
    outer_mass = []
    half_r=radius/2
    
    for i,j in zip(r_array,am):
        if i < half_r:
            enclosed_mass.append(j)
        else:
            outer_mass.append(j)
    enclosed_mass = np.array(enclosed_mass)
    outer_mass = np.array(outer_mass)
    mass=np.sum(enclosed_mass)
    mass_out=np.sum(outer_mass)
    print('mass in 1/2 radius = ',mass)
    print('mass out 1/2 radius = ',mass_out)
    print('num particles in 1/2 radius',enclosed_mass.shape[0])
    print('num particles out 1/2 radius',outer_mass.shape[0])

    r_orig_array=np.argsort(r_array)
    ascending_r = np.sort(r_array)

    half_mass=np.sum(am)/2
    #print(f'mass * 1/2 = {half_mass}')
    half_mass_array = np.array([])
    for x,y in zip(ascending_r,r_orig_array):
        if np.sum(half_mass_array) < half_mass:
            half_mass_array = np.append(half_mass_array,am[y])

    print('half mass = ', np.sum(half_mass_array))
    print('num particles in inner half mass ', half_mass_array.shape[0])
    
    GAM = 5/3
    print(f"GAMMA = {GAM}")

    # writes all data to bestfit.sph
    with open("bestfit.sph", "w") as f:
        anoteqfrac = 0

        # divides all values by the mass of the bin to convert back to the final values
        for index in range(nbinsbf):
            ammrhoavg[index] = ammrhoavg[index] / amavg[index]
            uiavg[index] = uiavg[index] / amavg[index]
            rhoavg[index] = rhoavg[index] / amavg[index]
            ravg[index] = ravg[index] / amavg[index]
            jrotavg[index] = np.linalg.norm(jrotavg[index]) / amavg[index]
            tavg[index] = tavg[index] / amavg[index]
            pavg[index] = pavg[index] / amavg[index]

            # apressure = uiavg[index] * rhoavg[index] * (GAM - 1)

            # updates composition data:
            h1avg[index] /= amavgcomp[index]
            he3avg[index] /= amavgcomp[index]
            he4avg[index] /= amavgcomp[index]
            c12avg[index] /= amavgcomp[index]
            n14avg[index] /= amavgcomp[index]
            o16avg[index] /= amavgcomp[index]
            ne20avg[index] /= amavgcomp[index]
            mg24avg[index] /= amavgcomp[index]

            # writes all values to bestfit.sph. If there is a NaN value
            # anywhere in the code, all values for that bin will be skipped
            if (not np.isnan(ammrhoavg[index] * amasstot * munit) or not np.isnan(ravg[index] *\
                runit) or not np.isnan(pavg[index] * gravconst * (munit / runit ** 2) ** 2) or not\
                np.isnan(rhoavg[index] * munit / runit ** 3) or not np.isnan(uiavg[index] * \
                gravconst * munit / runit) or not np.isnan(np.linalg.norm(jrotavg[index]) * \
                runit**2 / (np.sqrt(runit**3/(gravconst * munit))))) and \
                (pavg[index] != 0 or uiavg[index] != 0 or tavg[index] != 0):
                f.write('{:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e}\n'.format(
                ammrhoavg[index] * amasstot * munit,
                ravg[index] * runit,
                # apressure * gravconst * (munit / runit ** 2) ** 2,
                pavg[index],
                rhoavg[index] * munit / runit ** 3,
                uiavg[index] * gravconst * munit / runit,
                np.linalg.norm(jrotavg[index]) * runit**2 / tunit,
                tavg[index],
                h1avg[index], he3avg[index], he4avg[index], c12avg[index], n14avg[index],
                o16avg[index], ne20avg[index], mg24avg[index]))
            else:
                pass

            if index != 0:
                if uiavg[index]/rhoavg[index]**(GAM-1) < uiavg[index - 1]/rhoavg[index]**(GAM-1) \
                    and anoteqfrac == 0:  #This finds where A=(GAM-1)*u/rho^(GAM-1) stops increasing outward
                    anoteqfrac = 1 - 0.5 * (ammrhoavg[index]+ammrhoavg[index-1])
                    print(f"FRACTION NOT IN EQUILIBRIUM ={anoteqfrac} ")

        #apressure = aimax * rhomin * (GAM - 1)

        # writes the data for the outermost layers of the star
        f.write("{:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e}\n".format(amasstot * munit, radius * runit,
                                                       # apressure * gravconst * (munit/runit**2)**2,
                                                       pmax,
                                                       rhomin * munit/runit**3,
                                                       aimax * gravconst * munit/runit,
                                                       np.linalg.norm(jrotmax) * runit**2 / tunit,
                                                       tmax,
                                                       h1max, he3max, he4max, c12max, n14max, o16max, ne20max, mg24max
                                                       ))

        f.close()

    with open('sph_star.dat','w') as f:
        f.write(f'Mass:                             {np.sum(am)}\n')
        f.write(f'Radius:                           {radius}\n')
        f.write(f'Particles:                        {ntot}\n')
        f.write(f'Particles for inner 50% of mass:  {half_mass_array.shape[0]}\n')

        f.close()

    print("Leaving BESTFIT")
