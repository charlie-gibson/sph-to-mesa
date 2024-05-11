"""
This code attempts to to find the number of neighbors by calling a kdtree
This information can then be used to help smooth out the profiles that are being created
and can help to remove any jagged lines in the profiles
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

# wendland c4 kernel
def kernel(r,h):
    u = r/h
    sigW = 495. / (256. * 3.14159265358979323846264338327950288419716939937510582097494459237)
    q = u/2.
    result = sigW * (1. - q)**6 * (35./3. * q**2 + 6.*q + 1.)

    return result

# function that finds the neighbors of every particle
def neighbors(readit_data,comp_data,sph_input,component_val,component_data,neos=1):

<<<<<<< HEAD
    from bound_particles import bound_particle_data

    from get_temperature import get_temperature
    from eos_func import read_eos
    from eos_func import useeostable


    # constants used for calculations
    gravconst = 6.67390e-08 # m^3 / (kg s^2)
    munit = 1.9891e33 # kg
    runit = 6.9599e10 # m
    tunit = np.sqrt(runit**3/(gravconst * munit)) # s
    boltz = 1.380658e-16 # erg/kelvin
    crad = 2.997924580e10 # cm/sec  NOTE: crad has a different meaning in MESA
    planck = 6.6260755e-27 # gram cm^2/sec
    crad2 = crad**2 # cm^2/s^2
    sigma = np.pi**2*boltz*(boltz*2*np.pi/planck)**3/60/crad2 #cgs
    arad = 4.0*sigma/crad #cgs
    qconst = 1.5*boltz/arad #cgs
=======
    gravconst = 6.67390e-08 # m^3 / (kg s^2)
    munit = 1.9891e33 # kg
    runit = 6.9599e10 # m
    tunit = np.sqrt(runit**3/(gravconst * munit))
>>>>>>> 80b51c09b6d14c8c80b29d98c4f2b29bcd6c776d

    # gets the necessary values from the output file
    x=readit_data['x']
    y=readit_data['y']
    z=readit_data['z']
    vx=readit_data['vx']
    vy=readit_data['vy']
    vz=readit_data['vz']
    hp=readit_data['hp']
    rho=readit_data['rho']
    u=readit_data['u']
    am=readit_data['am']
    meanmolecular=readit_data['meanmolecular']
    n=readit_data['ntot']
    h1=np.array(comp_data[0])
    he3=np.array(comp_data[1])
    he4=np.array(comp_data[2])
    c12=np.array(comp_data[3])
    n14=np.array(comp_data[4])
    o16=np.array(comp_data[5])
    ne20=np.array(comp_data[6])
    mg24=np.array(comp_data[7])

    jrot=np.zeros(n)
    temp=np.zeros(n)
    press=np.zeros(n)

    if component_val != 0:
        breadit_data,bcomposition,j1,j2,j3=bound_particle_data(readit_data,component_data,comp_data,jrot,temp,press,component_val)

        bx=breadit_data['x']
        by=breadit_data['y']
        bz=breadit_data['z']
        bvx=breadit_data['vx']
        bvy=breadit_data['vy']
        bvz=breadit_data['vz']
        bhp=breadit_data['hp']
        brho=breadit_data['rho']
        bu=breadit_data['u']
        bam=breadit_data['am']
        bmeanmolecular=breadit_data['meanmolecular']
        bn=breadit_data['ntot']
        bh1=bcomposition[0]
        bhe3=bcomposition[1]
        bhe4=bcomposition[2]
        bc12=bcomposition[3]
        bn14=bcomposition[4]
        bo16=bcomposition[5]
        bne20=bcomposition[6]
        bmg24=bcomposition[7]

    else:

        bx=readit_data['x']
        by=readit_data['y']
        bz=readit_data['z']
        bvx=readit_data['vx']
        bvy=readit_data['vy']
        bvz=readit_data['vz']
        bhp=readit_data['hp']
        brho=readit_data['rho']
        bu=readit_data['u']
        bam=readit_data['am']
        bmeanmolecular=readit_data['meanmolecular']
        bn=readit_data['ntot']
        bh1=composition[0]
        bhe3=composition[1]
        bhe4=composition[2]
        bc12=composition[3]
        bn14=composition[4]
        bo16=composition[5]
        bne20=composition[6]
        bmg24=composition[7]

    inds = breadit_data['indices']
    mtot = np.sum(bam)

    try:
        ind = np.argmax(brho)
    except:
        print('No Particles Bound to Star: Stopping Code')
        raise SystemExit
    true_ind = inds[ind]

    # redefines the origin
    xc = bx[ind]
    yc = by[ind]
    zc = bz[ind]

    print('com = ',xc,yc,zc)

    # calculates velocity of the origin
    vxc=bvx[ind]
    vyc=bvy[ind]
    vzc=bvz[ind]

    print('v_com = ',vxc,vyc,vzc)

    # udpates all positions to have the correct radial distance from the origin

    xr=np.array(x-xc)
    yr=np.array(y-yc)
    zr=np.array(z-zc)

    ind = true_ind

    # double checks that the origin is at (0, 0, 0)

    print(f'r: {xr[ind]}, {yr[ind]}, {zr[ind]}')

        # redefines velocities accounting for origin drift
    vxr=np.array(vx-vxc)
    vyr=np.array(vy-vyc)
    vzr=np.array(vz-vzc)

    # double checks that the velocity is at (0, 0, 0)
    print(f'v: {vxr[ind]}, {vyr[ind]}, {vzr[ind]}')

    # creates three dimensional radius and velocity vectors
    rvec=np.zeros((n,3))
    vvec=np.zeros((n,3))

    for j in range(n):
        rvec[j,0]=xr[j]
        rvec[j,1]=yr[j]
        rvec[j,2]=zr[j]

        vvec[j,0]=vxr[j]
        vvec[j,1]=vyr[j]
        vvec[j,2]=vzr[j]

    # double checks that the radius and velocity at the center of mass are still (0, 0, 0)
    print(f'r: {rvec[ind]}')
    print(f'v: {vvec[ind]}')
    rcrossvind = np.cross(rvec[ind], vvec[ind])

    # calculates specific angular momentum
    rcrossv = np.cross(rvec,vvec)
    print(r"r$\times$ v: "+f"{rcrossv}")

    # creates arrays for each component of angular momentum
    jx = np.zeros(n)
    jy = np.zeros(n)
    jz = np.zeros(n)

    for i in range(n):
        jx[i] = rcrossv[i][0]
        jy[i] = rcrossv[i][1]
        jz[i] = rcrossv[i][2]

    # creates new empty lists to contain smoothed values
    jrotnew=np.zeros(n)
    tempnew=np.zeros(n)
    pnew=np.zeros(n)

    vxnew=np.zeros(n)
    vynew=np.zeros(n)
    vznew=np.zeros(n)

    jxnew = np.zeros(n)
    jynew = np.zeros(n)
    jznew = np.zeros(n)

    
    print('CALCULATING NEIGHBORS')
    # creates a numpy array with shape (n_particles, 3) with the three being for each dimension
    neighbors=np.zeros((n,3))

    # places each particle's position value into the array
    for k in range(n):
        neighbors[k,0]=x[k]
        neighbors[k,1]=y[k]
        neighbors[k,2]=z[k]

    #print(my_array2)

    # uses scipy.cKDTree (c meaning the language) to find every particle's neighbors
    tree=cKDTree(neighbors)

    # indicies is a list of lists with each nested list containing all the (sorted) indices that
    # are in a smoothing length, hp, from the particle
    indices=[np.array(tree.query_ball_point(point,2*h,return_sorted=False)) for point,h in zip(neighbors,hp)]

    print('FINISHED CALCULATING NEIGHBORS\n')

    # updates composition for each particle
    h1new=np.zeros(n)
    he3new=np.zeros(n)
    he4new=np.zeros(n)
    c12new=np.zeros(n)
    n14new=np.zeros(n)
    o16new=np.zeros(n)
    ne20new=np.zeros(n)
    mg24new=np.zeros(n)

    #### slower ####
    # pindex=0 # particle number
    # for i in indices: # i is the array of all neighbors of the particle
    # # indices is a list of lists where i is just one particle's neighbors
    #     r=np.sqrt((x[i]-x[pindex])**2+(y[i]-y[pindex])**2+(z[i]-z[pindex])**2)
    #     h=hp[pindex]
    #     kernelval=am[i]*kernel(r,h)/(h**3*rho[pindex])

    #     h1new[pindex]+=np.sum(h1[i]*kernelval)
    #     he3new[pindex]+=np.sum(he3[i]*kernelval)
    #     he4new[pindex]+=np.sum(he4[i]*kernelval)
    #     c12new[pindex]+=np.sum(c12[i]*kernelval)
    #     n14new[pindex]+=np.sum(n14[i]*kernelval)
    #     o16new[pindex]+=np.sum(o16[i]*kernelval)
    #     ne20new[pindex]+=np.sum(ne20[i]*kernelval)
    #     mg24new[pindex]+=np.sum(mg24[i]*kernelval)

    #     jxnew[pindex]+=np.sum(jx[i]*kernelval)
    #     jynew[pindex]+=np.sum(jy[i]*kernelval)
    #     jznew[pindex]+=np.sum(jz[i]*kernelval)

    #     pindex+=1

    ### faster ####
    sum_h1 = np.zeros_like(h1new)
    sum_he3 = np.zeros_like(he3new)
    sum_he4 = np.zeros_like(he4new)
    sum_c12 = np.zeros_like(c12new)
    sum_n14 = np.zeros_like(n14new)
    sum_o16 = np.zeros_like(o16new)
    sum_ne20 = np.zeros_like(ne20new)
    sum_mg24 = np.zeros_like(mg24new)

    sum_jx = np.zeros_like(jxnew)
    sum_jy = np.zeros_like(jynew)
    sum_jz = np.zeros_like(jznew)

    # Iterate over particles and their neighbors
    for pindex, neighbors in enumerate(indices):
        # Calculate distances between particle pindex and its neighbors
        dx = x[neighbors] - x[pindex]
        dy = y[neighbors] - y[pindex]
        dz = z[neighbors] - z[pindex]
        distances = np.sqrt(dx**2 + dy**2 + dz**2)

        h = hp[pindex]
        kernelval = am[neighbors] * kernel(distances, h) / (h**3 * rho[pindex])

        # Update the sums
        sum_h1[pindex] = np.sum(h1[neighbors] * kernelval)
        sum_he3[pindex] = np.sum(he3[neighbors] * kernelval)
        sum_he4[pindex] = np.sum(he4[neighbors] * kernelval)
        sum_c12[pindex] = np.sum(c12[neighbors] * kernelval)
        sum_n14[pindex] = np.sum(n14[neighbors] * kernelval)
        sum_o16[pindex] = np.sum(o16[neighbors] * kernelval)
        sum_ne20[pindex] = np.sum(ne20[neighbors] * kernelval)
        sum_mg24[pindex] = np.sum(mg24[neighbors] * kernelval)

<<<<<<< HEAD
        sum_jx[pindex] = np.sum(jx[neighbors] * kernelval)
        sum_jy[pindex] = np.sum(jy[neighbors] * kernelval)
        sum_jz[pindex] = np.sum(jz[neighbors] * kernelval)

    # Update the new arrays
    h1new += sum_h1
    he3new += sum_he3
    he4new += sum_he4
    c12new += sum_c12
    n14new += sum_n14
    o16new += sum_o16
    ne20new += sum_ne20
    mg24new += sum_mg24

    jxnew += sum_jx
    jynew += sum_jy
    jznew += sum_jz

    rcrossv = np.cross(rvec, vvec)

    rcrossvind2 = np.cross(rvec[ind], vvec[ind])
    print(f'r: {rvec[ind]}')
    # double checks that the angular momentum at the center is still (0, 0, 0) after the smoothing process
    print(f'rcrossv ind: {rcrossvind}')
    print(f'rcrossv ind 2: {rcrossvind2}')

    for i in range(n):
        rcrossv[i,0] = jxnew[i]
        rcrossv[i,1] = jynew[i]
        rcrossv[i,2] = jznew[i]

    jrotnew = rcrossv
=======
    # updates composition for each particle
    # print(h1[0], h1[1], h1[2])
    h1new=np.zeros(n)
    he3new=np.zeros(n)
    he4new=np.zeros(n)
    c12new=np.zeros(n)
    n14new=np.zeros(n)
    o16new=np.zeros(n)
    ne20new=np.zeros(n)
    mg24new=np.zeros(n)

    pindex=0 # particle number
    for i in indices: # i is the array of all neighbors of the particle
    # indices is a list of lists where i is just one particle's neighbors
        r=np.sqrt((x[i]-x[pindex])**2+(y[i]-y[pindex])**2+(z[i]-z[pindex])**2)
        h=hp[pindex]
        kernelval=am[i]*kernel(r,h)/(h**3*rho[pindex])
        h1new[pindex]+=np.sum(h1[i]*kernelval)
        he3new[pindex]+=np.sum(he3[i]*kernelval)
        he4new[pindex]+=np.sum(he4[i]*kernelval)
        c12new[pindex]+=np.sum(c12[i]*kernelval)
        n14new[pindex]+=np.sum(n14[i]*kernelval)
        o16new[pindex]+=np.sum(o16[i]*kernelval)
        ne20new[pindex]+=np.sum(ne20[i]*kernelval)
        mg24new[pindex]+=np.sum(mg24[i]*kernelval)

#        jrotnew[pindex]+=np.sum(rcrossv[i]*kernelval)

#        tempnew[pindex]+=np.sum(kernelval*temp)

        pindex+=1

    #print(h1new[0],h1new[1],h1new[2])

    # ratioH1=h1new/h1
    # ratioHe3=he3new/he3
    # ratioHe4=he4new/he4
    # ratioC12=c12new/c12
    # ratioN14=n14new/n14
    # ratioO16=o16new/o16
    # ratioNe20=ne20new/ne20
    # ratioMg24=mg24new/mg24

    # f,ax=plt.subplots()
    # ax.scatter(range(n),ratioH1,label='H1',marker='.',s=0.5)
    # ax.scatter(range(n),ratioHe3,label='He3',marker='.',s=0.5)
    # ax.scatter(range(n),ratioHe4,label='He4',marker='.',s=0.5)
    # ax.scatter(range(n),ratioC12,label='C12',marker='.',s=0.5)
    # ax.scatter(range(n),ratioN14,label='N14',marker='.',s=0.5)
    # ax.scatter(range(n),ratioO16,label='C12',marker='.',s=0.5)
    # ax.scatter(range(n),ratioNe20,label='Ne20',marker='.',s=0.5)
    # ax.scatter(range(n),ratioMg24,label='Mg24',marker='.',s=0.5)
    # ax.set_xlabel('Particle')
    # ax.set_ylabel('new composition / old composition')
    # ax.set_yscale('log')
    # ax.legend()
    # plt.show()
>>>>>>> 80b51c09b6d14c8c80b29d98c4f2b29bcd6c776d
    
           
    composition=[h1new,he3new,he4new,c12new,n14new,o16new,ne20new,mg24new]

<<<<<<< HEAD
    return composition,jrotnew,tempnew,pnew
=======
    jrotnew,tempnew,pnew=smooth_temp(readit_data,neos,composition,indices,component_val,component_data,sph_input)

    #jrotnew=jrotnew*runit**2/tunit
    print(tempnew)

    return composition,jrotnew,tempnew,pnew

def smooth_temp(readit_data,neos,composition,indices,component_val,component_data,sph_input):
    from bound_particles import bound_particle_data

    from get_temperature import get_temperature
    from eos_func import read_eos
    from eos_func import useeostable

    print('SMOOTHING TEMPERATURE')

    # gets the necessary values from the output file
    x=readit_data['x']
    y=readit_data['y']
    z=readit_data['z']
    vx=readit_data['vx']
    vy=readit_data['vy']
    vz=readit_data['vz']
    hp=readit_data['hp']
    rho=readit_data['rho']
    u=readit_data['u']
    am=readit_data['am']
    meanmolecular=readit_data['meanmolecular']
    n=readit_data['ntot']
    h1=composition[0]
    he3=composition[1]
    he4=composition[2]
    c12=composition[3]
    n14=composition[4]
    o16=composition[5]
    ne20=composition[6]
    mg24=composition[7]

    jrot=np.zeros(n)
    temp=np.zeros(n)
    press=np.zeros(n)

    if component_val != 0:
        breadit_data,bcomposition,j1,j2,j3=bound_particle_data(readit_data,component_data,composition,jrot,temp,press,component_val)

        bx=breadit_data['x']
        by=breadit_data['y']
        bz=breadit_data['z']
        bvx=breadit_data['vx']
        bvy=breadit_data['vy']
        bvz=breadit_data['vz']
        bhp=breadit_data['hp']
        brho=breadit_data['rho']
        bu=breadit_data['u']
        bam=breadit_data['am']
        bmeanmolecular=breadit_data['meanmolecular']
        bn=breadit_data['ntot']
        bh1=bcomposition[0]
        bhe3=bcomposition[1]
        bhe4=bcomposition[2]
        bc12=bcomposition[3]
        bn14=bcomposition[4]
        bo16=bcomposition[5]
        bne20=bcomposition[6]
        bmg24=bcomposition[7]

    else:

        bx=readit_data['x']
        by=readit_data['y']
        bz=readit_data['z']
        bvx=readit_data['vx']
        bvy=readit_data['vy']
        bvz=readit_data['vz']
        bhp=readit_data['hp']
        brho=readit_data['rho']
        bu=readit_data['u']
        bam=readit_data['am']
        bmeanmolecular=readit_data['meanmolecular']
        bn=readit_data['ntot']
        bh1=composition[0]
        bhe3=composition[1]
        bhe4=composition[2]
        bc12=composition[3]
        bn14=composition[4]
        bo16=composition[5]
        bne20=composition[6]
        bmg24=composition[7]
        

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

    temp = np.zeros(n)
    pgas = np.zeros(n)
    prad = np.zeros(n)
    ptot = np.zeros(n)
    # number is ntot - 1 to account for one particle being a black hole
    # I would like to find a way to make this line more general to avoid having to edit
    # the code for each new run
    rhocgs = np.zeros(n)
    ucgs = np.zeros(n)
    qval = np.zeros(n)
    rval = np.zeros(n)

    if neos == 1:
        rhocgs = rho*munit/runit**3
        # print(rhocgs)
        for i in range(n):
            ucgs[i] = u[i]*gravconst*munit/runit
        # print(ucgs)
        qval = qconst*rhocgs/meanmolecular
        # print(qval)
        rval = -ucgs*rhocgs/arad
        # print(rval)

        for nval in range(n):
            try:
                temp[nval] = get_temperature(qval[nval], rval[nval])
                pgas[nval] = rhocgs[nval] * boltz * temp[nval] / (meanmolecular[nval])
                # prad[n] = arad * temp[n]**4/3
                ptot[nval] = pgas[nval] # + prad[n]
            except:
                print("ERROR CALCULATING TEMPERATURE OR PRESSURE")

    elif neos == 2:
        eosfile=sph_input['eosfile']
        eos_data = read_eos(eosfile)
        # print(eos_data)
        for i in range(n):
            try:
                rhocgs[i] = rho[i]*munit/runit**3
                ucgs[i] = u[i]*gravconst*munit/runit
                # qval[i] = qconst*rhocgs[i]/meanmolecular[i]
                # rval[i] = -ucgs[i]*rhocgs[i]/arad

                temp[i] = useeostable(eos_data=eos_data, ucgs=ucgs[i], rhocgs=rhocgs[i], xxx=h1[i], which=0)
                ptot[i] = useeostable(eos_data=eos_data, ucgs=ucgs[i], rhocgs=rhocgs[i], xxx=h1[i], which=2)
                pgas[i] = ptot[i] - arad * temp[i]**4/3
                ptot[i] = pgas[i]
            except:
                print("problem using eostable", i, rho[i]*munit/runit**3, u[i]*gravconst*munit/runit, useeostable(eos_data=eos_data, ucgs=ucgs[i], rhocgs=rhocgs[i], xxx=h1[i], which=0))
                raise SystemExit

    print(temp)
    breadit_data,bcomposition,j1,j2,j3=bound_particle_data(readit_data,component_data,composition,jrot,temp,press,component_val)
    print(j2)

    mtot=np.sum(bam)

    xc=np.sum(bx*bam)/mtot
    yc=np.sum(by*bam)/mtot
    zc=np.sum(bz*bam)/mtot

    print('c.o.m.=',xc,yc,zc)

    vxc=np.sum(bvx*bam)/mtot
    vyc=np.sum(bvy*bam)/mtot
    vzc=np.sum(bvz*bam)/mtot

    print('v_com=',vxc,vyc,vzc)

    xr=np.array(x-xc)
    yr=np.array(y-yc)
    zr=np.array(z-zc)

    vxr=np.array(vx-vxc)
    vyr=np.array(vy-vyc)
    vzr=np.array(vz-vzc)

    rvec=np.zeros((n,3))
    vvec=np.zeros((n,3))
    
    for j in range(n):
        rvec[j,0]=xr[j]
        rvec[j,1]=yr[j]
        rvec[j,2]=zr[j]

        vvec[j,0]=vxr[j]
        vvec[j,1]=vyr[j]
        vvec[j,2]=vzr[j]

    rcrossv=np.cross(rvec,vvec)

    jrotmagn=np.zeros(n)
    for k in range(n):
        jrotmagn[k]=np.linalg.norm(rcrossv[k])
    #jrotmagn=np.linalg.norm(rcrossv)

    jrotnew=np.zeros(n)
    tempnew=np.zeros(n)
    pnew=np.zeros(n)

    pindex=0
    for i in indices:
        r=np.sqrt((x[i]-x[pindex])**2+(y[i]-y[pindex])**2+(z[i]-z[pindex])**2)
        h=hp[pindex]
        kernelval=am[i]*kernel(r,h)/(h**3*rho[pindex])
        tempnew[pindex]+=np.sum(temp[i]*kernelval)
        jrotnew[pindex]+=np.sum(jrotmagn[i]*kernelval)
        pnew[pindex]+=np.sum(pgas[i]*kernelval)

        pindex+=1

    print(tempnew)

    return jrotnew,tempnew,pnew
>>>>>>> 80b51c09b6d14c8c80b29d98c4f2b29bcd6c776d
