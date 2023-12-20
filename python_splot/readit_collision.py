"""
Hans Vanderzyden 4/24/23 Python Splot Code for Reading in Binary out.sph Files in StarSmasher

This python code reads in StarSmasher data from out.sph files in the current directory. 
This collects the particle data for each particle in the simulation, with the last particle being a black hole core particle (if applicable).

This was created for a collision run, may work for relaxation runs <---------- FIGURE OUT

---------------------------------------------------------------------------------------------------
Update: Charles Gibson 05/20/2023
The functionality is the same as the original readit.py file. However, I made several quality of
life features.

First, I reworked the code to use the mmap module. This optimizes the readin of the file to
significantly speed up the process.

Second, I placed all the data into a dictionary called data. This dictionary is then returned, so
the data can be passed to other codes when called upon in the splot.py interface code.

Note that the original code was written by Hans Vanderzyden and ChatGPT was used to help with a
lot of the conversion
"""

def readit_collision(nnit, iform):
    
    # imports the necessary modules for the binary file readin
    import mmap # used for optimization
    import os
    import numpy as np
    import math

    # globalizes all variables so they can be accessed by other codes when this is imported elsewhere
    global x, y, z, am, hp, rho, vx, vy, vz, vxdot, \
    vydot, vzdot, u, udot, grpot, meanmolecular, cc, divv, \
    aa, bb, dd
    x = [] # x position of the particle
    y = [] # y position of the particle
    z = [] # z position of the particle
    am = [] # mass of the particle
    hp = [] # smoothing length of the particle
    rho = [] # density of the particle
    vx = [] # x velocity of the particle
    vy = [] # y velocity of the particle
    vz = [] # z velocity of the particle
    vxdot = [] # x acceleration of the particle
    vydot = [] # y acceleration of the particle
    vzdot = [] # z acceleration of the particle
    u = [] # specific internal energy
    udot = [] # d/dt(specific internal energy)
    grpot = [] # gravitational potential of the particle
    meanmolecular = [] # mean molecular weight of the particle
    cc = [] # integer labeling which parent star the particle came from
    divv = []
    aa = []
    bb = []
    dd = []

    # used to ensure correct naming of the file
    if nnit < 10:
        FNAME = f'out000{nnit}.sph'
    elif nnit < 100:
        FNAME = f'out00{nnit}.sph'
    elif nnit < 1000:
        FNAME = f'out0{nnit}.sph'
    else:
        FNAME = f'out{nnit}.sph'

    print(f'Looking for {FNAME}')

    fileexists = os.path.isfile(FNAME)

    # reads in the file(s) if it (they) exist
    if fileexists:

        # reads the file in as binary
        with open(FNAME, 'rb') as f:

            # uses mmap to read in the header data
            # these values correspond to the input data from sph.input in the StarSmasher run
            with mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ) as mm:
                offset = 4
                ntot, nnopt = np.frombuffer(mm, dtype=np.int32, count=2, offset=offset)
                offset += 2 * np.dtype(np.int32).itemsize
                hco, hfloor, sep0, tf, dtout = np.frombuffer(mm, dtype=np.float64, count=5, offset=offset)
                offset += 5 * np.dtype(np.float64).itemsize
                nout, nit = np.frombuffer(mm, dtype=np.int32, count=2, offset=offset)
                offset += 2 * np.dtype(np.int32).itemsize
                t = np.frombuffer(mm, dtype=np.float64, count=1, offset=offset)[0]
                offset += np.dtype(np.float64).itemsize
                nav = np.frombuffer(mm, dtype=np.int32, count=1, offset=offset)[0]
                offset += np.dtype(np.int32).itemsize
                alpha, beta, tjumpahead = np.frombuffer(mm, dtype=np.float64, count=3, offset=offset)
                offset += 3 * np.dtype(np.float64).itemsize
                ngr, nrelax = np.frombuffer(mm, dtype=np.int32, count=2, offset=offset)
                offset += 2 * np.dtype(np.int32).itemsize
                trelax, dt, omega2 = np.frombuffer(mm, dtype=np.float64, count=3, offset=offset)
                offset += 3 * np.dtype(np.float64).itemsize
                ncooling = np.frombuffer(mm, dtype=np.int32, count=1, offset=offset)[0]
                offset += np.dtype(np.int32).itemsize
                erad = np.frombuffer(mm, dtype=np.float64, count=1, offset=offset)[0]
                offset += np.dtype(np.float64).itemsize
                ndisplace = np.frombuffer(mm, dtype=np.int32, count=1, offset=offset)[0]
                offset += np.dtype(np.int32).itemsize
                displacex, displacey, displacez = np.frombuffer(mm, dtype=np.float64, count=3, offset=offset)
                offset += 3 * np.dtype(np.float64).itemsize
                junk = np.frombuffer(mm, dtype=np.float64, count=1, offset=offset)[0]
                print("HEADERS READ IN")

                # prints number of particles
                global NCHK
                NCHK = ntot
                print('NCHK=', NCHK)

                m = 0

                # assigns values for each particle
                for i in range(ntot):
                    offset += np.dtype(np.float64).itemsize
                    xvar, yvar, zvar, amvar, hpvar, rhovar, vxvar, vyvar, vzvar, vxdotvar,\
                        vydotvar, vzdotvar, uvar, udotvar, grpotvar, meanmolecularvar \
                        = np.frombuffer(mm, dtype=np.float64, count=16, offset=offset)
                    offset += np.dtype(np.float64).itemsize * 16
                    ccvar = np.frombuffer(mm, dtype=np.int32, count=1, offset=offset)[0]
                    offset += np.dtype(np.int32).itemsize
                    divvvar, aavar, bbvar, ddvar = np.frombuffer(mm, dtype=np.float64, count=4, offset=offset)
                    offset += np.dtype(np.float64).itemsize * 4
                    junk = np.frombuffer(mm, dtype=np.float64, count=1, offset=offset)[0]

                    # adds the particle data to each list intialized above
                    # check if i can speed this up without append
                    if rhovar != 0:
                        x.append(xvar)
                        y.append(yvar)
                        z.append(zvar)
                        am.append(amvar)
                        hp.append(hpvar)
                        rho.append(rhovar)
                        vx.append(vxvar)
                        vy.append(vyvar)
                        vz.append(vzvar)
                        vxdot.append(vxdotvar)
                        vydot.append(vydotvar)
                        vzdot.append(vzdotvar)
                        u.append(uvar)
                        udot.append(udotvar)
                        grpot.append(grpotvar)
                        meanmolecular.append(meanmolecularvar)
                        cc.append(ccvar)
                        divv.append(divvvar)
                        aa.append(aavar)
                        bb.append(bbvar)
                        dd.append(ddvar)
                    else:
                        m += 1

                print("DATA READ IN")

                # converts all data to numpy arrays
                x = np.array(x)
                y = np.array(y)
                z = np.array(z)
                am = np.array(am)
                hp = np.array(hp)
                rho = np.array(rho)
                vx = np.array(vx)
                vy = np.array(vy)
                vz = np.array(vz)
                vxdot = np.array(vxdot)
                vydot = np.array(vydot)
                vzdot = np.array(vzdot)
                u = np.array(u)
                udot = np.array(udot)
                grpot = np.array(grpot)
                meanmolecular = np.array(meanmolecular)
                cc = np.array(cc)
                divv = np.array(divv)
                aa = np.array(aa)
                bb = np.array(bb)
                dd = np.array(dd)

    # this runs if there is no more matching file requested
    else:
        print('Ran out of output files')
    
    ntot -= m

    # updates all velocities based on their accelerations
    # to corerct the leapfrog integration technique
    vx -= 0.5*dt*vxdot
    vy -= 0.5*dt*vydot
    vz -= 0.5*dt*vzdot
    u = [u[i] - 0.5*dt*udot[i] if u[i] != 0.0 else u[i] for i in range(ntot)]
    t -= 0.5*dt
    
    print("VELOCITIES UPDATED")

    # used for corotating binary: more work/understanding needed for translation
    if omega2 != 0.0 and nrelax == 2:
        theta = math.sqrt(omega2) * t
        print('omega2=', omega2, theta)
        print('CALCULATION WAS IN ROTATING FRAME: UNDOING ROTATION!')

        for i in range(int(ntot)):
            # rotation is counterclockwise
            xold = x[i]
            yold = y[i]
            x[i] = xold * np.cos(theta) - yold * np.sin(theta)
            y[i] = xold * np.sin(theta) + yold * np.cos(theta)
            vxold = vx[i]
            vyold = vy[i]
            vx[i] = vxold * np.cos(theta) - vyold * np.sin(theta) - np.sqrt(omega2) * x[i]
            vy[i] = vxold * np.sin(theta) + vyold * np.cos(theta) + np.sqrt(omega2) * y[i]

    try:
        with open('sph.start1u', 'rb') as f:
            f.seek(4) # skips the first 4 bytes of the file (there is some sort of junk data, idk why?)
            n1 = np.fromfile(f, dtype=np.int32, count=1, sep="")
            n1=int(n1)
            cc1val=cc[0]
            

        with open('sph.start2u','rb') as f:
            f.seek(4) # skips the first 4 bytes of the file (there is some sort of junk data, idk why?)
            n2 = np.fromfile(f, dtype=np.int32, count=1, sep="")
            n2=int(n2)
            cc2val=cc[-1]
            
        print(f'n1 = {n1}     n2 = {n2}')
            
        if ntot!=n1+n2:
            print('discrepancy between ntot and n1+n2')

    except:
        pass

    # places all lists of data into a dictionary to be passed to other files
    try:
        data = {
            'ntot':int(ntot),
            'x':x,
            'y':y,
            'z':z,
            'am':am,
            'hp':hp,
            'rho':rho,
            'vx':vx,
            'vy':vy,
            'vz':vz,
            'vxdot':vxdot,
            'vydot':vydot,
            'vzdot':vzdot,
            'u':u,
            'udot':udot,
            'grpot':grpot,
            'meanmolecular':meanmolecular,
            'cc':cc,
            'divv':divv,
            'aa':aa,
            'bb':bb,
            'dd':dd,
            'dtout':dtout,
            'dt':dt,
            'n1':n1,
            'n2':n2,
            'cc1val':cc1val,
            'cc2val':cc2val
        }
    except:
        data = {
            'ntot':int(ntot),
            'x':x,
            'y':y,
            'z':z,
            'am':am,
            'hp':hp,
            'rho':rho,
            'vx':vx,
            'vy':vy,
            'vz':vz,
            'vxdot':vxdot,
            'vydot':vydot,
            'vzdot':vzdot,
            'u':u,
            'udot':udot,
            'grpot':grpot,
            'meanmolecular':meanmolecular,
            'cc':cc,
            'divv':divv,
            'aa':aa,
            'bb':bb,
            'dd':dd,
            'dtout':dtout,
            'dt':dt
        }

    return data
    
