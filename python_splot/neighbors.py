"""
This code attempts to to find the number of neighbors by calling a kdtree
This information can then be used to help smooth out the profiles that are being created
and can help to remove any jagged lines in the profiles
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

# def kernel(r,h):
#     constant=1/(np.pi*h**3)
#     result=np.zeros(len(r))

#     condition1=r>2*h
#     condition2=r<=h
#     condition3=(r>h) & (r<=2*h)

#     #print(condition1.shape, condition2.shape, condition3.shape)
#     #print(condition1, condition2, condition3)

#     result=np.where(condition1,0,result)
#     result=np.where(condition2,1/(np.pi*h**3)*(1-3/2*(r/h)**2+3/4*(r/h)**3),result)
#     result=np.where(condition3,1/(np.pi*h**3)*(1/4*(2-(r/h))**3),result)
#     #result[~condition1 & ~condition2] = (np.pi*h**3)[~condition1 & ~condition2] * (1/4*(2-r[~condition1 & ~condition2]/h[~condition1 & ~condition2])**3)
    
#     return result

# wendland c4 kernel
def kernel(r,h):
    u = r/h
    sigW = 495. / (256. * 3.14159265358979323846264338327950288419716939937510582097494459237)
    q = u/2.
    result = sigW * (1. - q)**6 * (35./3. * q**2 + 6.*q + 1.)

    return result

# function that finds the neighbors of every particle
def neighbors(readit_data,comp_data,sph_input,component_val,component_data,neos=1):

    gravconst = 6.67390e-08 # m^3 / (kg s^2)
    munit = 1.9891e33 # kg
    runit = 6.9599e10 # m
    tunit = np.sqrt(runit**3/(gravconst * munit))

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
    #indices=[np.array(tree.query(point,h)) for point,h in zip(neighbors,hp)]
    #print(indices)
    #print(indices[0])
    # for i in indices:
    #     i.pop()

    #print(indices)
    
    #raise SystemExit

    print('FINISHED CALCULATING NEIGHBORS\n')

    # # checks the neighbor calculation by comparing density
    # rhonew=np.zeros(n)
    
    # pindex=0
    # for i in indices:
    #     # pindex is the index of the particle as what corresponds to the original arrays

    #     # the following declarations allow for numpy vectorization of the density calculation

    #     # i is a numpy array of the indices that are neighbors of particle i
    #     #i=np.array(i)
    #     # finds the distance between the particle and its neighbors
    #     r=np.sqrt((x[i]-x[pindex])**2+(y[i]-y[pindex])**2+(z[i]-z[pindex])**2)
    #     # finds the mass of the neighboring particles
    #     # finds the smoothing length of the neighboring particles
    #     h=hp[pindex]

    #     # debugging and output
    #     #if pindex%10000==0:
    #         #print(r, type(r),r.shape)
    #         #print(m, type(m),m.shape)
    #         #print(h, type(h),h.shape)
    #         #print(r)
    #         #print(pindex)

    #     # the density of the particle associated with pindex is the sum of the mass of each
    #     # neighboring particle times the kernel value based on the distance from the particle
    #     # and the particle's smoothing length
    #     rhonew[pindex]+=np.sum(am[i]*kernel(r,h))/h**3
    #     pindex+=1


    #    # for k in i:
    #         #r=np.sqrt((x[i]-x[k])**2+(y[i]-y[k])**2+(z[i]-z[k])**2)
    #     #    r=np.sqrt((x[k]-x[pindex])**2+(y[k]-y[pindex])**2+(z[k]-z[pindex])**2)
    #         #print(r,hp[pindex])
    #      #   rhonew += m[k] * kernel(r,hp[pindex])

    # #print(rho)
    # #print(rhonew)

    # ratio=rhonew/rho
    
    # f,ax=plt.subplots()
    # ax.scatter(range(n),ratio,s=0.5,marker='.',color='black')
    # ax.set_xlabel('Particle')
    # ax.set_ylabel(r'$\frac{\rho_{new}}{\rho_{old}}$')
    # plt.show()

    # # raise SystemExit

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
    
    # raise SystemExit
    # for i in indices:
    #     pindex=indices.index(i)
    #     for k in i:
    #         r=np.sqrt((x[k]-x[pindex])**2+(y[d]-y[pindex])**2+(z[k]-z[pindex])**2)
    #         print(r)
    #         h1new[i]+=h1[i]*m[k]*kernel(r,h[i])/rho[i]
    #         he3new[i]+=he3[i]*m[k]*kernel(r,h[i])/rho[i]
    #         he4new[i]+=he4[i]*m[k]*kernel(r,h[i])/rho[i]
    #         c12new[i]+=c12[i]*m[k]*kernel(r,h[i])/rho[i]
    #         n14new[i]+=n14[i]*m[k]*kernel(r,h[i])/rho[i]
    #         o16new[i]+=o16[i]*m[k]*kernel(r,h[i])/rho[i]
    #         ne20new[i]+=ne20[i]*m[k]*kernel(r,h[i])/rho[i]
    #         mg24new[i]+=mg24[i]*m[k]*kernel(r,h[i])/rho[i]
           
    composition=[h1new,he3new,he4new,c12new,n14new,o16new,ne20new,mg24new]

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
