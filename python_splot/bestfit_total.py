"""
This script calls the py file from splot_directories.

Using the read in data, this code sorts every particle based on the density
into many bins. It then calculates certain values such as density, pressure,
specific internal energy, mass, radius, angular momentum, and composition for
each bin. Finally, it bins these into a lower resolution set of bins and writes
them to a new file called bestfit.sph.

The original Fortran file states:
Analyze the results of a collision run and output important
numbers for a summary table

Charles Gibson
Allegheny College
Department of Physics
05/26/2023
"""

# data is read in using readit.py and passed through splot.py
def bestfit(data, comp_data):

    from readit import readit
    import numpy as np

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
    divv = data['divv']
    aa = data['aa']
    bb = data['bb']
    dd = data['dd']
    dt = data['dt'] * 0.018445 * 24 * 60 * 60 # to get dt in seconds
    
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

    
    # h1 = np.array(h1).astype(float)
    # he3 = np.array(he3).astype(float)
    # he4 = np.array(he4).astype(float)
    # c12 = np.array(c12).astype(float)
    # n14 = np.array(n14).astype(float)
    # o16 = np.array(o16).astype(float)
    # ne20 = np.array(ne20).astype(float)
    # mg24 = np.array(mg24).astype(float)

    print("DATA PASSED TO FUNCTION")
    
    # variables
    nbins = 1601 # number of bins for sorting data
    nbinsbf = 100 # number of bins for the best fit data: higher values means higher resolution
    alogrhomin = -7.9 # 10^-10 is the lowest density
    alogrhomax = 2.1 # 10^10 is the highest density

    # constants used for calculations
    gravconst = 6.67390e-08 # m^3 / (kg s^2)
    munit = 1.9891e33 # kg
    runit = 6.9599e10 # m
    tunit = 0.018445 * 24 * 60 * 60 # days * hr * min * sec

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
    aiavg = np.zeros(nbinsbf)
    ravg = np.zeros(nbinsbf)
    amavg = np.zeros(nbinsbf)

    # loops through every particle
    for i in range(ntot):
        alogrhoi = np.log10(rho[i]) # finds the log base 10 of each density value
        # determines which bin to put the particle
        nbini = int((alogrhoi-alogrhomin) * (nbins-1.0)/(alogrhomax - alogrhomin) + 1.0)

        # checks if the bin is outside of the range of possible bins
        if(nbini > nbins or nbini <  1):
            print("PROBLEM 3")

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
    ammrhoavg = np.zeros(nbinsbf+1)
    amavg = np.zeros(nbinsbf+1)
    aiavg = np.zeros(nbinsbf+1)
    rhoavg = np.zeros(nbinsbf+1)
    ravg = np.zeros(nbinsbf+1)
    jrotavg = np.zeros((nbinsbf+1, 3))
    radius = 0

    # composition arrays
    h1avg = np.zeros(nbinsbf+1)
    he3avg = np.zeros(nbinsbf+1)
    he4avg = np.zeros(nbinsbf+1)
    c12avg = np.zeros(nbinsbf+1)
    n14avg = np.zeros(nbinsbf+1)
    o16avg = np.zeros(nbinsbf+1)
    ne20avg = np.zeros(nbinsbf+1)
    mg24avg = np.zeros(nbinsbf+1)

    # loops over every particle
    for i in range(ntot):
        alogrhoi = np.log10(rho[i]) # finds the log of the density of each particle
        nbini = (alogrhoi - alogrhomin) * (nbins - 1) / (alogrhomax - alogrhomin) + 1 # index of the bin that the particle is in
        amonmrho = amrho[int(nbini)] / amasstot # mass fraction of the particle based on the bin it falls into
        index = int(amonmrho * nbinsbf + 0.99999)
        if index > nbinsbf:
            print("Warning: index > nbinsbf", index)
            index = nbinsbf
        
        if index < 1:
            print("Warning: index < 1", index)
            index = 1

        # note that all of these multiply by mass but are divided by another mass later on
        amavg[index] += am[i] # adds the mass for each particle to its corresponding index
        aiavg[index] += u[i] * am[i] # this gives the internal energy for each index
        rhoavg[index] += rho[i] * am[i] # this gives the density times mass for each index
        rvec = np.array([x[i] - xc, y[i] - yc, z[i] - zc]) # this provides the r vector of the particle
        vvec = np.array([(vx[i]-vxc), (vy[i]-vyc), (vz[i]-vzc)]) # this gives the velocity vector of the particle
        r = np.linalg.norm(rvec) # this is the length of the r vector from the center of mass
        v = np.linalg.norm(vvec) # this is the speed of the v vector relative to the center of mass
        ravg[index] += r * am[i] # this gives the distance times the mass of the particle
        jrotavg[index] += am[i] * np.cross(rvec, vvec) # this gives the angular momentum times the mass of the particle

        # composition data
        h1avg[index] += h1[i] * am[i]
        he3avg[index] += he3[i] * am[i]
        he4avg[index] += he4[i] * am[i]
        c12avg[index] += c12[i] * am[i]
        n14avg[index] += n14[i] * am[i]
        o16avg[index] += o16[i] * am[i]
        ne20avg[index] += ne20[i] * am[i]
        mg24avg[index] += mg24[i] * am[i]

        # finds these values for the outermost particle
        if r > radius:
            radius = r
            aimax = u[i]
            rhomin = rho[i]
            jrotmax = am[i] * np.cross(rvec, vvec)

            h1max = h1[i]
            he3max = he3[i]
            he4max = he4[i]
            c12max = c12[i]
            n14max = n14[i]
            o16max = o16[i]
            ne20max = ne20[i]
            mg24max = mg24[i]

        ammrhoavg[index] = ammrhoavg[index] + amonmrho * am[i]

    GAM = 5/3
    print(f"GAMMA = {GAM}")

    # writes all data to bestfit.sph
    with open("bestfit.sph", "w") as f:
        anoteqfrac = 0

        # divides all values by the mass of the bin to convert back to the final values
        for index in range(nbinsbf):
            ammrhoavg[index] = ammrhoavg[index] / amavg[index]
            aiavg[index] = aiavg[index] / amavg[index]
            rhoavg[index] = rhoavg[index] / amavg[index]
            ravg[index] = ravg[index] / amavg[index]
            jrotavg[index] = np.linalg.norm(jrotavg[index]) / amavg[index]

            apressure = aiavg[index] * rhoavg[index] * (GAM - 1)

            # updates composition data:
            h1avg[index] /= amavg[index]
            he3avg[index] /= amavg[index]
            he4avg[index] /= amavg[index]
            c12avg[index] /= amavg[index]
            n14avg[index] /= amavg[index]
            o16avg[index] /= amavg[index]
            ne20avg[index] /= amavg[index]
            mg24avg[index] /= amavg[index]

            # writes all values to bestfit.sph. If there is a NaN value
            # anywhere in the code, it will be skipped
            if not np.isnan(ammrhoavg[index] * amasstot * munit) or not np.isnan(ravg[index] *\
                runit) or not np.isnan(apressure * gravconst * (munit / runit ** 2) ** 2) or not\
                np.isnan(rhoavg[index] * munit / runit ** 3) or not np.isnan(aiavg[index] * \
                gravconst * munit / runit) or not np.isnan(np.linalg.norm(jrotavg[index]) * \
                runit**2 / (np.sqrt(runit**3/(gravconst * munit)))):
                f.write('{:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e}\n'.format(
                ammrhoavg[index] * amasstot * munit,
                ravg[index] * runit,
                apressure * gravconst * (munit / runit ** 2) ** 2,
                rhoavg[index] * munit / runit ** 3,
                aiavg[index] * gravconst * munit / runit,
                np.linalg.norm(jrotavg[index]) * runit**2 / (np.sqrt(runit**3/(gravconst * munit))),
                h1avg[index], he3avg[index], he4avg[index], c12avg[index], n14avg[index],
                o16avg[index], ne20avg[index], mg24avg[index]))
            else:
                pass

            if index != 1:
                if aiavg[index] < aiavg[index - 1] and anoteqfrac == 0:
                    anoteqfrac = 1 - 0.5 * (ammrhoavg[index]+ammrhoavg[index-1])
                    print(f"FRACTION NOT IN EQUILIBRIUM ={anoteqfrac} ")

        apressure = aimax * rhomin * (GAM - 1)

        # writes the data for the outermost layers of the star
        f.write("{:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e}\n".format(amasstot * munit, radius * runit,
                                                       apressure * gravconst * (munit/runit**2)**2,
                                                       rhomin * munit/runit**3,
                                                       aimax * gravconst * munit/runit,
                                                       np.linalg.norm(jrotmax) * runit**2 / (np.sqrt(runit**3/(gravconst * munit))),
                                                       h1max, he3max, he4max, c12max, n14max, o16max, ne20max, mg24max
                                                       ))

        f.close()

    print("Leaving BESTFIT")
