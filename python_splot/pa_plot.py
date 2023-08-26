import numpy as np
import matplotlib.pyplot as plt
import math
from get_temperature import get_temperature

# outfile = int(input("out file number: "))

# data = readit_collision(outfile,4)

def pa_plot(data, profile_num, outfile):

    # saves necessary data from readit in splot.py
    ntot = int(data['ntot']) # number of particles
    x = data['x'] # x position
    y = data['y'] # y position
    z = data['z'] # z position
    am = data['am'] # mass
    hp = data['hp'] # smoothing length
    rho = data['rho'] # density
    u = data['u'] # specific internal energy
    meanmolecular = data['meanmolecular'] # mean molecular weight

    # constants are used to calculate the values in cgs
    gravconst = 6.67390e-08 # m^3 / (kg s^2)
    munit = 1.9891e33 # kg
    runit = 6.9599e10 # m
    tunit = np.sqrt(runit**3/(gravconst*munit))
    boltz = 1.380658e-16 # erg/kelvin
    crad = 2.997924580e10 # cm/sec  NOTE: crad has a different meaning in MESA
    planck = 6.6260755e-27 # gram cm^2/sec
    crad2 = crad**2 # cm^2/s^2
    sigma = np.pi**2*boltz*(boltz*2*np.pi/planck)**3/60/crad2 #cgs
    arad = 4.0*sigma/crad #cgs
    qconst = 1.5*boltz/arad #cgs

    # calculates the temperature using get_temperature() function used by SPH
    temp = np.zeros(ntot)
    pgas = np.zeros(ntot)
    prad = np.zeros(ntot)
    ptot = np.zeros(ntot)
    rhocgs = np.zeros(ntot)
    ucgs = np.zeros(ntot)
    qval = np.zeros(ntot)
    rval = np.zeros(ntot)

    rhocgs = rho*munit/runit**3
    for i in range(ntot):
        ucgs[i] = u[i]*gravconst*munit/runit
    qval = qconst*rhocgs/meanmolecular
    rval = -ucgs*rhocgs/arad

    for n in range(ntot):
        temp[n] = get_temperature(qval[n], rval[n])
        pgas[n] = rhocgs[n] * boltz * temp[n] / (meanmolecular[n])
        prad[n] = arad * temp[n]**4/3
        ptot[n] = pgas[n] + prad[n]

    for j in ptot:
        j=j/((munit*runit/tunit**2)/(munit**2))

    # determines the radius and A (P/rho^(5/3)) values
    r = np.sqrt(x**2 + y**2 + z**2)
    A = ptot/rhocgs**(5/3)

    # reads in the profile used in the SPH relaxation
    profilename = profile_num

    profilemass = []
    profileradius = []
    profilepressure = []
    profiledensity = []
    profileA = []

    with open(profilename,'r') as f:
        n=0
        # skips the first several lines which don't contain information about pressure, density
        for line in f:
            if n < 6:
                n+=1
            else:
                # appends the values to the new lists. Change the indices if the formatting of the profile is different
                vals = line.split()
                profileradius.append(10**float(vals[2]))
                profiledensity.append(10**float(vals[4]))
                profilepressure.append(10**float(vals[5]))

    # calculates A
    for i in range(len(profiledensity)):
        profileA.append(profilepressure[i]/profiledensity[i]**(5/3))

    # reads in the data from the column file (col*.sph) to get the neighbor number and accelerations
    neighbornumber = []
    radialg = []
    ahydro = []

    # rounds the column to the next highest multiple of 10
    colname = int(math.floor(outfile/10)*10)

    # same naming convention as in readit_collision.py
    if colname < 10:
        FNAME = f'col000{colname}.sph'
    elif colname < 100:
        FNAME = f'col00{colname}.sph'
    elif colname < 1000:
        FNAME = f'col0{colname}.sph'
    else:
        FNAME = f'col{colname}.sph'

    # appends the neighbor number, radial g, and a_hydro to the three lists
    with open(FNAME,'r') as f:
        for line in f:
            line=line.split()
            neighbornumber.append(int(line[7]))
            radialg.append(float(line[9]))
            ahydro.append(float(line[8]))

    # plots all data in a subplot similar to that of the supermongo format
    f,axs=plt.subplots(4,2,gridspec_kw={'height_ratios':[1,1,1,1]})
    f.subplots_adjust(hspace=0,wspace=0.5)

    axs[0,0].scatter(r,A,c='black',s=0.25,marker='.')
    axs[0,0].plot(profileradius,profileA,c='red',linewidth=0.75)
    axs[0,0].set_yscale('log')
    axs[0,0].set_ylabel(r'$P/{\rho^{5/3}}$',fontsize=7)
    axs[0,0].set_xlim(0,np.max(r))
    axs[0,0].set_ylim(np.min(A),np.max(A))

    axs[1,0].scatter(r,ptot,c='black',s=0.25,marker='.')
    axs[1,0].plot(profileradius,profilepressure,c='red',linewidth=0.75)
    axs[1,0].set_yscale('log')
    axs[1,0].set_ylabel(r'${P}$',fontsize=7)
    axs[1,0].set_xlim(0,np.max(r))
    axs[1,0].set_ylim(np.min(ptot),np.max(ptot))

    axs[2,0].scatter(r,rhocgs,c='black',s=0.25,marker='.')
    axs[2,0].plot(profileradius,profiledensity,c='red',linewidth=0.75)
    axs[2,0].set_yscale('log')
    axs[2,0].set_ylabel(r'${\rho}$',fontsize=7)
    axs[2,0].set_xlim(0,np.max(r))
    axs[2,0].set_ylim(np.min(rhocgs),np.max(rhocgs))

    axs[0,1].scatter(r,am,c='black',s=0.25,marker='.')
    axs[0,1].set_ylabel(r'${m_i}$',fontsize=7)
    axs[0,1].set_yscale('log')

    axs[1,1].scatter(r,hp,c='black',s=0.25,marker='.')
    axs[1,1].set_ylabel(r'${h_i}$',fontsize=7)

    axs[2,1].scatter(r,neighbornumber,c='black',s=0.25,marker='.')
    axs[2,1].set_ylabel(r'${N_N}$',fontsize=7)

    axs[3,1].scatter(r,radialg,c='blue',s=0.25,marker='.')
    axs[3,1].scatter(r,ahydro,c='lime',s=0.25,marker='.')
    axs[3,1].set_ylabel(r'${radial g}$, $a_{hydro}$',fontsize=7)

    f.delaxes(axs[3,0])

    plt.show()
