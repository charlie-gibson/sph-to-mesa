"""
This script calls the py file from splot_directories.

This code determines gravitational bound components.

Charles Gibson and James Lombardi
Allegheny College
Department of Physics
"""

# data is read in using readit.py and passed through splot.py
def compbest3(nout, data):
    global firstt, icomp

    from readit_collision import readit_collision
    import numpy as np
    from calccom import calccom

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

    print("DATA RETRIEVED")

    #enth = u
    enth = np.zeros(ntot)

    icomp = np.empty(ntot)
    icomp = [1] * ntot

    #try:
    #    print(firstt)
    #except NameError:
    #    print("This the first time in this routine.")
    #    firstt = False
    #    ccparents = [cc[0]]  # Assign the first element of cc to ccparents[0]
    #    numparents = 1
    #    for i in range(ntot):
    #        j = 0
    #        while j < numparents:
    #            if cc[i] == ccparents[j]:
    #                break
    #            j += 1
    #        if j == numparents:
    #            numparents += 1
    #            ccparents.append(cc[i])
    #        icomp[i] = j
    #else:
    #    print("This isn't the first time in compbest3.")
    #print(icomp)

    ncomp1 = 0
    ncomp2 = 0
    ncomp3 = 0
    ncomp4 = 0
    rhomax1 = 0.0
    rhomax2 = 0.0
    rhomax3 = 0.0
    rhomax4 = 0.0
    irhomax1 = 0
    irhomax2 = 0
    irhomax3 = 0
    irhomax4 = 0

    for i in range(ntot):
        if icomp[i] == 1:
            if rho[i] > rhomax1:
                rhomax1 = rho[i]
                irhomax1 = i
            ncomp1 += 1
        elif icomp[i] == 2:
            if rho[i] > rhomax2:
                rhomax2 = rho[i]
                irhomax2 = i
            ncomp2 += 1
        elif icomp[i] == 3:
            if rho[i] > rhomax3:
                rhomax3 = rho[i]
                irhomax3 = i
            ncomp3 += 1
        elif icomp[i] == 4:
            if rho[i] > rhomax4:
                rhomax4 = rho[i]
                irhomax4 = i
            ncomp4 += 1

    print("irhomaxes =", irhomax1, irhomax2, irhomax3, irhomax4)
    print("number of component 1 particles in initial guess =", ncomp1)
    print("number of component 2 particles in initial guess =", ncomp2)
    print("number of component 3 particles in initial guess =", ncomp3)
    print("number of component 4 particles in initial guess =", ncomp4)

    am1,r1,v1,am2,r2,v2,am3,r3,v3,am4 = calccom(ntot,am,x,y,z,vx,vy,vz,icomp)
    
    #print(am1,r1,v1,am2,r2,v2,am3,r3,v3)
    print('nit nchng m1 m2 m3')
    print(0,ncomp1+ncomp2+ncomp3,am1,am2,am3)

    # Iterate over component determination
    ncit = 0
    NCITMAX = 100
    nchng = 1
    swapcase = 0
    while nchng != 0 and ncit < NCITMAX:
        nchng = 0
        for i in range(ntot):
            ri = np.array([x[i], y[i], z[i]])
            vi = np.array([vx[i], vy[i], vz[i]])
            # Calculate "total enthalpy" BI1 to component 1
            if am1 > am[i]:
                ri1 = np.linalg.norm(ri - r1)
                v2i1 = np.linalg.norm(vi - v1) ** 2
                bi1 = 0.5 * v2i1 + enth[i] - am1 / ri1
            else:
                bi1 = 1e30
                
            # calculate "total enthalpy" bi2 to component 2
            if am2 > am[i]:
                ri2 = np.linalg.norm(ri - r2)
                v2i2 = np.linalg.norm(vi - v2) ** 2
                bi2 = 0.5 * v2i2 + enth[i] - am2 / ri2
            else:
                bi2 = 1e30
                 
            # calculate "total enthalpy" bi3 to component 3
            if am3 > am[i]:
                ri3 = np.linalg.norm(ri - r3)
                v2i3 = np.linalg.norm(vi - v3) ** 2
                bi3 = 0.5 * v2i3 + enth[i] - am3 / ri3
            else:
                bi3 = 1e30
                     
            # assign particle to one of 3 components
            if bi1 < 0 and (bi2 >= 0 or ri1 < ri2) and (bi3 >= 0 or ri1 < ri3):
                ic = 1  # particle belongs to component 1
            elif bi2 < 0 and (bi1 >= 0 or ri2 < ri1) and (bi3 >= 0 or ri2 < ri3):
                ic = 2  # particle belongs to component 2
            elif bi3 < 0 and (bi1 >= 0 or ri3 < ri1) and (bi2 >= 0 or ri3 < ri2):
                ic = 3  # particle belongs to component 3
            else:
                ic = 4  # particle is unbound (component 4)
                        
            if ic != icomp[i]:
                # particle assignment has changed
                # adjust components
                icomp[i] = ic
                nchng += 1
        am1,r1,v1,am2,r2,v2,am3,r3,v3,am4 = calccom(ntot,am,x,y,z,vx,vy,vz,icomp)
        ncit += 1
        
        print(ncit,nchng,am1,am2,am3)

        if ncit >= NCITMAX:
           print('COMPBEST3: NO CONVERGENCE ???')
           raise SystemExit

    if nout <= 9999:
        fname = f'comp{nout:04d}.sph'
    else:
        fname = f'comp{nout:05d}.sph'
    print(f"Writing to file {fname}")
    with open(fname, 'w') as file:
        file.write('x      y      z      component\n')
        for i in range(ntot):
            file.write(f"{x[i]} {y[i]} {z[i]} {icomp[i]}\n")

    print("Leaving COMPBEST3")
