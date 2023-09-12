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
def neighbors(readit_data,comp_data):

    # gets the necessary values from the output file
    x=readit_data['x']
    y=readit_data['y']
    z=readit_data['z']
    hp=readit_data['hp']
    rho=readit_data['rho']
    am=readit_data['am']
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
    print(h1[0], h1[1], h1[2])
    h1new=np.zeros(n)
    he3new=np.zeros(n)
    he4new=np.zeros(n)
    c12new=np.zeros(n)
    n14new=np.zeros(n)
    o16new=np.zeros(n)
    ne20new=np.zeros(n)
    mg24new=np.zeros(n)

    pindex=0
    for i in indices:
        r=np.sqrt((x[i]-x[pindex])**2+(y[i]-y[pindex])**2+(z[i]-z[pindex])**2)
        h=hp[pindex]
        # h1i=np.array(h1[i])
        # he3i=np.array(he3[i])
        # he4i=np.array(he4[i])
        # c12i=np.array(c12[i])
        # n14i=np.array(n14[i])
        # o16i=np.array(o16[i])
        # ne20i=np.array(ne20[i])
        # mg24=np.array(mg24[i])
        # if pindex==0:
        #     print(pindex, h1[i], am[i], kernel(r,h))
        h1new[pindex]+=np.sum(h1[i]*am[i]*kernel(r,h))/(rho[pindex]*h**3)
        he3new[pindex]+=np.sum(he3[i]*am[i]*kernel(r,h))/(rho[pindex]*h**3)
        he4new[pindex]+=np.sum(he4[i]*am[i]*kernel(r,h))/(rho[pindex]*h**3)
        c12new[pindex]+=np.sum(c12[i]*am[i]*kernel(r,h))/(rho[pindex]*h**3)
        n14new[pindex]+=np.sum(n14[i]*am[i]*kernel(r,h))/(rho[pindex]*h**3)
        o16new[pindex]+=np.sum(o16[i]*am[i]*kernel(r,h))/(rho[pindex]*h**3)
        ne20new[pindex]+=np.sum(ne20[i]*am[i]*kernel(r,h))/(rho[pindex]*h**3)
        mg24new[pindex]+=np.sum(mg24[i]*am[i]*kernel(r,h))/(rho[pindex]*h**3)
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

    return composition
