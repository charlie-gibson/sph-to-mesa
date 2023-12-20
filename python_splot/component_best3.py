"""
This script calls the py file from splot_directories.

This code determines gravitational bound components.

Charles Gibson and James Lombardi
Allegheny College
Department of Physics
"""
# data is read in using readit.py and passed through splot.py
def compbest3(nout,data,icomp,write=True):
    global firstt

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
    try:
        n1 = data['n1']
        n2 = data['n2']
    except:
        pass
    
    print("DATA RETRIEVED")

    #enth = u
    enth = np.zeros(ntot)

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

    ri=np.column_stack([x,y,z])
    vi=np.column_stack([vx,vy,vz])

    orig_icomp=icomp
    orig_ic1=ncomp1
    orig_ic2=ncomp2
    orig_ic3=ncomp3

    #print(ri,vi)

    # Iterate over component determination
    ncit = 0
    NCITMAX = 200
    nchng = 1
    swapcase = 0
    retries=0
    while nchng != 0 and ncit < NCITMAX:
        nchng = 0

        valid_am1 = am1 > am
        valid_am2 = am2 > am
        valid_am3 = am3 > am

        #print(valid_am1, valid_am2, valid_am3)

        ri1 = np.linalg.norm(ri-r1,axis=1)
        ri2 = np.linalg.norm(ri-r2,axis=1)
        ri3 = np.linalg.norm(ri-r3,axis=1)

        #print(f'ri1 = {ri1}')

        v2i1 = np.linalg.norm(vi-v1,axis=1)**2
        v2i2 = np.linalg.norm(vi-v2,axis=1)**2
        v2i3 = np.linalg.norm(vi-v3,axis=1)**2

        #print(0.5 * v2i1)
        #print(am1/ri1, am2/ri2, am3/ri3)

        bi1 = np.where(valid_am1, 0.5 * v2i1 + enth - am1 / ri1, 1e30)
        bi2 = np.where(valid_am2, 0.5 * v2i2 + enth - am2 / ri2, 1e30)
        bi3 = np.where(valid_am3, 0.5 * v2i3 + enth - am3 / ri3, 1e30)

        #print(bi1,bi2,bi3)

        valid_bi1 = bi1 < 0
        valid_bi2 = bi2 < 0
        valid_bi3 = bi3 < 0

        #print(np.where(valid_bi1))

        valid_condition_1 = valid_bi1 & ((bi2 >= 0) | (ri1 < ri2)) & ((bi3 >= 0) | (ri1 < ri3))
        valid_condition_2 = valid_bi2 & ((bi1 >= 0) | (ri2 < ri1)) & ((bi3 >= 0) | (ri2 < ri3))
        valid_condition_3 = valid_bi3 & ((bi1 >= 0) | (ri3 < ri1)) & ((bi2 >= 0) | (ri3 < ri2))

#        ic = np.where(valid_condition_1, 1, np.where(valid_condition_2, 2, np.where(valid_condition_3, 3, 4)))

        #print(ic)

#        nchng = np.count_nonzero(ic != icomp)  # Count the number of changed assignments
#        icomp = ic  # Update the icomp array with the new assignments

        ic = np.where(valid_condition_1, 1, np.where(valid_condition_2, 2, np.where(valid_condition_3, 3, 4)))

        nchng = np.count_nonzero(ic != icomp)  # Count the number of changed assignments
        icomp = ic  # Update the icomp array with the new assignments
        #print(nchng)

        am1,r1,v1,am2,r2,v2,am3,r3,v3,am4 = calccom(ntot,am,x,y,z,vx,vy,vz,icomp)
        ncit += 1
        
        print(ncit,nchng,am1,am2,am3)

        # merges the two component 1 and component 2 assignments if the star doesn't converge
        convergence=True
        if ncit >= NCITMAX and retries == 0:
           #print('COMPBEST3: NO CONVERGENCE ???')
           icomp=np.array(orig_icomp)
           if orig_ic1>orig_ic2 and orig_ic2<3000:
               icomp2=np.where(icomp==2)[0]
               #print(icomp2)
               #raise SystemExit
               icomp[icomp2]=1
               am1,r1,v1,am2,r2,v2,am3,r3,v3,am4 = calccom(ntot,am,x,y,z,vx,vy,vz,icomp)
           elif orig_ic2>orig_ic1 and orig_ic1<3000:
               icomp1=np.where(icomp==1)[0]
               #print(icomp1)
               icomp[icomp1]=2
               #print(icomp[np.where(icomp==1)[0]])
               am1,r1,v1,am2,r2,v2,am3,r3,v3,am4 = calccom(ntot,am,x,y,z,vx,vy,vz,icomp)
           else:
               print('TOO MANY PARTICLES TO MERGE STARS. EXITING CODE')
               raise SystemExit
           print('NO CONVERGENCE: MERGING STARS AND RESTARTING\n')
           retries+=1
           ncit=0
           #raise SystemExit
        elif ncit>=NCITMAX:
            print('COMPBEST3: NO CONVERGENCE')
            convergence=False

    if nout <= 9999:
        fname = f'comp{nout:04d}.sph'
    else:
        fname = f'comp{nout:05d}.sph'

    if write:
        print(f"Writing to file {fname}")
        with open(fname, 'w') as file:
            file.write('x      y      z      component\n')
            for i in range(ntot):
                file.write(f"{x[i]} {y[i]} {z[i]} {int(icomp[i])}\n")

    print("Leaving COMPBEST3")

    return icomp,convergence
