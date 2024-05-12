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

# reads in the collision
def readit_collision(nnit, iform=4, ascii = False):
    
    """
    nnit: output number
    iform: from original Fortran splot - not used - set to 4
    ascii: True if the data being read in is to be written to an ascii output file
    """

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
        m = 0
        with open(FNAME, 'rb') as f:

            # creates a memory map of the out*.sph file
            with mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ) as mm:
                # accounts for junk data at the beginning of the file
                offset = 4
        
                # Read header data: integer values have 32 bytes and float values have 64 bytes
                header_dtype = np.dtype([
                    ('ntot', np.int32), ('nnopt', np.int32), ('hco', np.float64), ('hfloor', np.float64),
                    ('sep0', np.float64), ('tf', np.float64), ('dtout', np.float64), ('nout', np.int32),
                    ('nit', np.int32), ('t', np.float64), ('nav', np.int32), ('alpha', np.float64),
                    ('beta', np.float64), ('tjumpahead', np.float64), ('ngr', np.int32), ('nrelax', np.int32),
                    ('trelax', np.float64), ('dt', np.float64), ('omega2', np.float64), ('ncooling', np.int32),
                    ('erad', np.float64), ('ndisplace', np.int32), ('displacex', np.float64), ('displacey', np.float64),
                    ('displacez', np.float64), ('junk', np.float64)
                ])
                header_data = np.frombuffer(mm, dtype=header_dtype, count=1, offset=offset)

                # Extracts values from header data
                ntot = header_data['ntot'][0]
                nnopt = header_data['nnopt'][0]
                hco = header_data['hco'][0]
                hfloor = header_data['hfloor'][0]
                sep0 = header_data['sep0'][0]
                tf = header_data['tf'][0]
                dtout = header_data['dtout'][0]
                nout = header_data['nout'][0]
                nit = header_data['nit'][0]
                t = header_data['t'][0]
                nav = header_data['nav'][0]
                alpha = header_data['alpha'][0]
                beta = header_data['beta'][0]
                tjumpahead = header_data['tjumpahead'][0]
                ngr = header_data['ngr'][0]
                nrelax = header_data['nrelax'][0]
                trelax = header_data['trelax'][0]
                dt = header_data['dt'][0]
                omega2 = header_data['omega2'][0]
                ncooling = header_data['ncooling'][0]
                erad = header_data['erad'][0]
                ndisplace = header_data['ndisplace'][0]
                displacex = header_data['displacex'][0]
                displacey = header_data['displacey'][0]
                displacez = header_data['displacez'][0]

                # creates a list of header data for later access
                # removes reliance on memory map
                header_data = [ntot,nnopt,hco,hfloor,sep0,tf,dtout,nout,nit,t,nav,alpha,\
                beta,tjumpahead,ngr,nrelax,trelax,dt,omega2,ncooling,erad,ndisplace,\
                displacex,displacey,displacez]

                # Batch read particle data
                particle_dtype = np.dtype([
                    ('x', np.float64), ('y', np.float64), ('z', np.float64), ('am', np.float64),
                    ('hp', np.float64), ('rho', np.float64), ('vx', np.float64), ('vy', np.float64),
                    ('vz', np.float64), ('vxdot', np.float64), ('vydot', np.float64), ('vzdot', np.float64),
                    ('u', np.float64), ('udot', np.float64), ('grpot', np.float64), ('meanmolecular', np.float64),
                    ('cc', np.int32), ('divv', np.float64), ('aa', np.float64), ('bb', np.float64),
                    ('dd', np.float64), ('junk', np.float64)
                ])
                # avoids having to loop over the entire file
                particle_data = np.frombuffer(mm, dtype=particle_dtype, count=ntot, offset=offset + np.dtype(header_dtype).itemsize)
                
                # Extract values from particle data
                x = particle_data['x']
                y = particle_data['y']
                z = particle_data['z']
                am = particle_data['am']
                hp = particle_data['hp']
                rho = particle_data['rho']
                vx = particle_data['vx']
                vy = particle_data['vy']
                vz = particle_data['vz']
                vxdot = particle_data['vxdot']
                vydot = particle_data['vydot']
                vzdot = particle_data['vzdot']
                u = particle_data['u']
                udot = particle_data['udot']
                grpot = particle_data['grpot']
                meanmolecular = particle_data['meanmolecular']
                cc = particle_data['cc']
                divv = particle_data['divv']
                aa = particle_data['aa']
                bb = particle_data['bb']
                dd = particle_data['dd']
                
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
    
    # ntot -= m

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

    # sees if there are initial stars in the collision to get the number of particles
    # from star 1 and star 2
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
<<<<<<< HEAD

=======
<<<<<<< HEAD

=======
            
>>>>>>> 80b51c09b6d14c8c80b29d98c4f2b29bcd6c776d
>>>>>>> 706dca6185b661b857fdf6d3337257d9d53c74ff
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
    # if there are not separate n1 and n2 values (relaxation, TDE, etc.)
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

    # returns header data for ascii output (splot option 0)
    if ascii:
        return data, header_data
    # otherwise, skips header_data being returned
    else:
        return data
