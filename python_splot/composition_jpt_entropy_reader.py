"""
This code analyzes the file bestfit2.sph to generate
a profile pair of density and specific thermal energy.

These values can be used to generate a mesa run using
the relax_composition_j_entropy test suite.

As of now, the composition will be used from the original
mesa run until we can get the composition in the output
from StarSmasher.

Charles Gibson
Allegheny College
Department of Physics
"""

def entropy_reader(mode):

    import numpy as np
    from scipy.interpolate import CubicSpline
    import matplotlib.pyplot as plt

    # reads in the data from bestfit.sph
    with open("bestfit.sph") as f:
        # headers
        # header_names = [ "mass", "radius", "pressure", "density", "specific thermal energy" "specific angular momentum"]

        # instantiates empty lists of mass and thermodynamic variables
        mass = []
        radius = []
        pressure = []
        density = []
        specThermEnergy = []
        jrotlist = []
        templist = []
        h1list = []
        he3list = []
        he4list = []
        c12list = []
        n14list = []
        o16list = []
        ne20list = []
        mg24list = []

        # adds each value of the output to the above lists
        for line in f:
            # splits the lines into each variable and appends them to the appropriate list
            try:
                values = line.split()
                mass.append(float(values[0].replace("D", "e")))
                radius.append(float(values[1].replace("D", "e")))
                pressure.append(float(values[2].replace("D", "e")))
                density.append(float(values[3].replace("D", "e")))
                specThermEnergy.append(float(values[4].replace("D", "e")))
                jrotlist.append(float(values[5].replace("D", "e")))
                templist.append(float(values[6].replace("D", "e")))
                h1list.append(float(values[7].replace("D", "e")))
                he3list.append(float(values[8].replace("D", "e")))
                he4list.append(float(values[9].replace("D", "e")))
                c12list.append(float(values[10].replace("D", "e")))
                n14list.append(float(values[11].replace("D", "e")))
                o16list.append(float(values[12].replace("D", "e")))
                ne20list.append(float(values[13].replace("D", "e")))
                mg24list.append(float(values[14].replace("D", "e")))
            except:
                pass

    # temporary until we can figure out how num_zones is calculated in MESA
    with open("/home/kce5466/sph-to-mesa/python_splot/angular_momentum.dat") as f:
        
        xq = []
        q=[]

        n = 0

        for line in f:
            if n == 0:
                num_zones = line.split()
                n += 1
            else:
                newxq = line.split()[0]
                xq.append(float(newxq))
                q.append(1-float(newxq))

    totalmass = 0
    xq_mass_fraction = []
    q_mass_fraction = []

    # the total mass should be the final line's mass
    # which we find with this line

    totalmass = float(mass[-1])
    # print("entropy reader says normalized totalmass =", totalmass)
    # print(len(q), len(mass))

    index = 0
    # calculates the total mass outside each "shell"
    for n in mass:
        massenclosed = n # each value in the mass list is how much mass is found from the center to that radius
        massoutside = totalmass - massenclosed
        massfraction = massoutside / totalmass
        # print(massfraction)
        xq_mass_fraction.append(massfraction) # mesa uses this version of the mass fraction
        q_mass_fraction.append(1-massfraction)
        # print(xq_mass_fraction[index])
        index += 1

    # calculates P/rho^(5/3)
    poverrho53 = []
    for i, j in zip(density, pressure):
        poverrho53.append(j/i**(5/3))

    # Spline can only interpolate with increasing values of q
    xq_mass_fraction.reverse()
    q_mass_fraction.reverse()
    radius.reverse()
    density.reverse()
    specThermEnergy.reverse()
    pressure.reverse()
    poverrho53.reverse()
    jrotlist.reverse()
    templist.reverse()
    h1list.reverse()
    he3list.reverse()
    he4list.reverse()
    c12list.reverse()
    n14list.reverse()
    o16list.reverse()
    ne20list.reverse()
    mg24list.reverse()

    # increases the resolution during low resolution areas in the star

#    resolved = False
#    while resolved == False:
        
#        for i in range(len(q_mass_fraction)):
#            try:
#                if q_mass_fraction[i+1]-q_mass_fraction[i]>=0.04:
#                    q_mass_fraction.insert(i+1,q_mass_fraction[i]+0.01)
#                    mass.insert(i+1,(mass[i]+mass[i+1])/2)
#                    radius.insert(i+1,(radius[i]+radius[i+1])/2)
#                    pressure.insert(i+1,(pressure[i]+pressure[i+1])/2)
#                    density.insert(i+1,(density[i]+density[i+1])/2)
#                    specThermEnergy.insert(i+1,(specThermEnergy[i]+specThermEnergy[i+1])/2)
#                    jrotlist.insert(i+1,(jrotlist[i]+jrotlist[i+1])/2)
#                    templist.insert(i+1,(templist[i]+templist[i+1])/2)
#                    h1list.insert(i+1,(h1list[i]+h1list[i+1])/2)
#                    he3list.insert(i+1,(he3list[i]+he3list[i+1])/2)
#                    he4list.insert(i+1,(he4list[i]+he4list[i+1])/2)
#                    c12list.insert(i+1,(c12list[i]+c12list[i+1])/2)
#                    n14list.insert(i+1,(n14list[i]+n14list[i+1])/2)
#                    o16list.insert(i+1,(o16list[i]+o16list[i+1])/2)
#                    ne20list.insert(i+1,(ne20list[i]+ne20list[i+1])/2)
#                    mg24list.insert(i+1,(mg24list[i]+mg24list[i+1])/2)
#                    poverrho53.insert(i+1,(poverrho53[i]+poverrho53[i+1])/2)
#                    print('updated xq list')
#                    break
#                else:
#                    resolved=True
#            except:
#                pass

#    print(q_mass_fraction)
    
    # adds the value of the first zone to a point below zero to ensure good boundary condition behavior

    xq_mass_fraction.append((xq_mass_fraction[-1]+1)/2)
    q_mass_fraction.append((q_mass_fraction[-1])/2)
    mass.append(mass[-1])
    radius.append(radius[-1])
    pressure.append(pressure[-1])
    density.append(density[-1])
    specThermEnergy.append(specThermEnergy[-1])
    jrotlist.append(jrotlist[-1])
    templist.append(templist[-1])
    h1list.append(h1list[-1])
    he3list.append(he3list[-1])
    he4list.append(he4list[-1])
    c12list.append(c12list[-1])
    n14list.append(n14list[-1])
    o16list.append(o16list[-1])
    ne20list.append(ne20list[-1])
    mg24list.append(mg24list[-1])
    poverrho53.append(poverrho53[-1])

    
    xq_mass_fraction.append(1.0)
    q_mass_fraction.append(0.0)
    mass.append(mass[-1])
    radius.append(radius[-1])
    pressure.append(pressure[-1])
    density.append(density[-1])
    specThermEnergy.append(specThermEnergy[-1])
    jrotlist.append(jrotlist[-1])
    templist.append(templist[-1])
    h1list.append(h1list[-1])
    he3list.append(he3list[-1])
    he4list.append(he4list[-1])
    c12list.append(c12list[-1])
    n14list.append(n14list[-1])
    o16list.append(o16list[-1])
    ne20list.append(ne20list[-1])
    mg24list.append(mg24list[-1])
    poverrho53.append(poverrho53[-1])
    
    xq_mass_fraction.append(1.01)
    q_mass_fraction.append(-0.01)
    mass.append(mass[-1])
    radius.append(radius[-1])
    pressure.append(pressure[-1])
    density.append(density[-1])
    specThermEnergy.append(specThermEnergy[-1])
    jrotlist.append(jrotlist[-1])
    templist.append(templist[-1])
    h1list.append(h1list[-1])
    he3list.append(he3list[-1])
    he4list.append(he4list[-1])
    c12list.append(c12list[-1])
    n14list.append(n14list[-1])
    o16list.append(o16list[-1])
    ne20list.append(ne20list[-1])
    mg24list.append(mg24list[-1])
    poverrho53.append(poverrho53[-1])

    xq_mass_fraction.append(1.02)
    q_mass_fraction.append(-0.02)
    mass.append(mass[-1])
    radius.append(radius[-1])
    pressure.append(pressure[-1])
    density.append(density[-1])
    specThermEnergy.append(specThermEnergy[-1])
    jrotlist.append(jrotlist[-1])
    templist.append(templist[-1])
    h1list.append(h1list[-1])
    he3list.append(he3list[-1])
    he4list.append(he4list[-1])
    c12list.append(c12list[-1])
    n14list.append(n14list[-1])
    o16list.append(o16list[-1])
    ne20list.append(ne20list[-1])
    mg24list.append(mg24list[-1])
    poverrho53.append(poverrho53[-1])

    # creates an equation for the data of mass fraction vs density,
    # mass fraction vs specific thermal energy, and mass fraction and angular momentum
    rho = CubicSpline(xq_mass_fraction, density)
    e = CubicSpline(xq_mass_fraction, specThermEnergy)
    p = CubicSpline(xq_mass_fraction, pressure)
    prho53 = CubicSpline(xq_mass_fraction, poverrho53)
    jrot = CubicSpline(xq_mass_fraction, jrotlist)
    temp = CubicSpline(xq_mass_fraction, templist)
    r = CubicSpline(xq_mass_fraction,radius)

    # uses linear interpolation for composition
    h1 = np.interp(x=xq,xp=xq_mass_fraction, fp=h1list)
    he3 = np.interp(x=xq,xp=xq_mass_fraction, fp=he3list)
    he4 = np.interp(x=xq,xp=xq_mass_fraction, fp=he4list)
    c12 = np.interp(x=xq,xp=xq_mass_fraction, fp=c12list)
    n14 = np.interp(x=xq,xp=xq_mass_fraction, fp=n14list)
    o16 = np.interp(x=xq,xp=xq_mass_fraction, fp=o16list)
    ne20 = np.interp(x=xq,xp=xq_mass_fraction, fp=ne20list)
    mg24 = np.interp(x=xq,xp=xq_mass_fraction, fp=mg24list)
    rcomp = np.interp(x=xq,xp=xq_mass_fraction, fp=radius)

    # now that spline is made, we need to write to the file
    xq_mass_fraction.reverse()
    q_mass_fraction.reverse()
    radius.reverse()
    density.reverse()
    specThermEnergy.reverse()
    pressure.reverse()
    poverrho53.reverse()
    jrotlist.reverse()
    templist.reverse()
    h1list.reverse()
    he3list.reverse()
    he4list.reverse()
    c12list.reverse()
    n14list.reverse()
    o16list.reverse()
    ne20list.reverse()
    mg24list.reverse()


    # creates the best fit equations to be graphed later
    densityfit = rho(xq)
    rfit = r(xq)
    specThermEnergyfit = e(xq)
    poverrho53fit = prho53(xq)
    pressurefit = p(xq)
    jrotfit = jrot(xq)
    tempfit = temp(xq)
    #h1fit = h1(xq)
    #he3fit = he3(xq)
    #he4fit = he4(xq)
    #c12fit = c12(xq)
    #n14fit = n14(xq)
    #o16fit = o16(xq)
    #ne20fit = ne20(xq)
    #mg24fit = mg24(xq)

    interp_data = {
        'xq':xq,
        'q':q,
        'r':rfit,
        'rho':densityfit,
        'P':pressurefit,
        'A':poverrho53fit,
        'T':tempfit,
        'jrot':jrotfit,
        'H1':h1,
        'He3':he3,
        'He4':he4,
        'C12':c12,
        'N14':n14,
        'O16':o16,
        'Ne20':ne20,
        'Mg24':mg24
    }
    
    mycompositionfile = open("composition.dat", "w")

    mycompositionfile.write(f"{num_zones[0]}    {8}\n")

    # evaluates all negative values of composition as 0

    
    
#    positiveh1 = [max(0, h1[k]) for k in q]
#    positiveh1 = [min(1, positiveh1[k]) for k in range(len(q))]
#    positivehe3 = [max(0, he3[k]) for k in q]
#    positivehe3 = [min(1, positivehe3[k]) for k in range(len(q))]
#    positivehe4 = [max(0, he4[k]) for k in q]
#    positivehe4 = [min(1, positivehe4[k]) for k in range(len(q))]
#    positivec12 = [max(0, c12[k]) for k in q]
#    positivec12 = [min(1, positivec12[k]) for k in range(len(q))]
#    positiven14 = [max(0, n14[k]) for k in q]
#    positiven14 = [min(1, positiven14[k]) for k in range(len(q))]
#    positiveo16 = [max(0, o16[k]) for k in q]
#    positiveo16 = [min(1, positiveo16[k]) for k in range(len(q))]
#    positivene20 = [max(0, ne20[k]) for k in q]
#    positivene20 = [min(1, positivene20[k]) for k in range(len(q))]
#    positivemg24 = [max(0, mg24[k]) for k in q]
#    positivemg24 = [min(1, positivemg24[k]) for k in range(len(q))]
        
    for k in range(len(xq)):
        mycompositionfile.write(f"{xq[k]}     {h1[k]}     {he3[k]}     {he4[k]}     {c12[k]}     {n14[k]}     {o16[k]}     {ne20[k]}     {mg24[k]}\n")
        

    # this writes the entropy.dat file to be used in mesa

    myentropyfile = open("entropy.dat", "w")

    # includes the number of zones at the top of the file
    myentropyfile.write(f"{num_zones[0]}\n")

    if mode == 'DT':

        # this format is for a density, temperature pair
        for j in xq:
            myentropyfile.write(f"{j}     {abs(rho(j))}     {abs(temp(j))}\n")

        myentropyfile.close()
    
    elif mode == 'PT':

        # print('length of q =', len(q))
        for j in xq:
            myentropyfile.write(f"{j}     {abs(p(j))}     {abs(temp(j))}\n")
        
        myentropyfile.close

    elif mode == 'DE':

        # this format is for a density, specific thermal energy pair
        # writes mass fraction, density, and specific thermal energy on one line per shell
        for j in xq:
            myentropyfile.write(f"{j}     {abs(rho(j))}     {abs(e(j))}\n")

        myentropyfile.close

    elif mode == 'DP':

        for j in xq:
            myentropyfile.write(f"{j}     {abs(rho(j))}     {abs(p(j))}\n")
        
        myentropyfile.close()

    # this writes the angular_momentum.dat file to be used in mesa

    myjrotfile = open("angular_momentum.dat", "w")

    # includes the number of zones at the top of the file
    myjrotfile.write(f"{num_zones[0]}\n")

    # writes mass fraction, density, and specific thermal energy on one line per shell
    for j in xq:
        myjrotfile.write(f"{j}     {jrot(j)}\n")

    myjrotfile.close()

    # calculates the total X, Y, and Z

    weightedX=0
    weightedY=0
    for i in range(len(xq)-1):
        xq_frac=xq[i+1]-xq[i]
        weightedX+=h1[i]*xq_frac
        weightedY+=he4[i]*xq_frac

    Z=1-weightedX-weightedY
    print(f'X: {weightedX}\nY: {weightedY}\nZ: {Z}')

    with open('sph_star.dat','a') as f:
        f.write(f'X:                                {weightedX}\n')
        f.write(f'Y:                                {weightedY}\n')
        f.write(f'Z:                                {Z}\n')
    
    # subplots the data for visual representation as a function of xq

    f,ax = plt.subplots(3,2,sharex=True,figsize=(8,8))

    f.subplots_adjust(hspace=0,wspace=0.5)
    
    ax[0,0].plot(xq,specThermEnergyfit,color='black')
    ax[0,0].set_ylabel(r'$u$')
    #ax[0,0].invert_xaxis()

    ax[1,0].plot(xq,densityfit,color='black')
    ax[1,0].set_ylabel(r'$\rho$')
    ax[1,0].set_yscale('log')
    #ax[1,0].invert_xaxis()

    ax[2,0].plot(xq,tempfit,color='black')
    ax[2,0].set_ylabel(r'$T$')
    ax[2,0].set_xlabel(r'$xq$')
    ax[2,0].set_xlim(1.1,-0.1)

    ax[0,1].plot(xq,poverrho53fit,color='black')
    #ax[0,1].scatter(q_mass_fraction,poverrho53,marker='x',color='black')
    ax[0,1].set_ylabel(r'${P/{\rho}^(5/3)}$')
    ax[0,1].set_yscale('log')
    #ax[0,1].invert_xaxis()
    
    ax[1,1].plot(xq,pressurefit,color='black')
    ax[1,1].set_ylabel(r'$P$')
    ax[1,1].set_yscale('log')
    #ax[1,1].invert_xaxis()

    ax[2,1].plot(xq,jrotfit,color='black')
    ax[2,1].set_ylabel(r'$j$')
    ax[2,1].set_xlabel(r'$xq$')
    ax[2,1].set_yscale('log')
    ax[2,1].set_xlim(1.1,-0.1)

    plt.show()

    # plots as a function or radius
    f,ax = plt.subplots(3,2,sharex=True,figsize=(8,8))

    f.subplots_adjust(hspace=0,wspace=0.5)
    
    ax[0,0].plot(rfit,specThermEnergyfit,color='black')
    ax[0,0].set_ylabel(r'$u$')
    #ax[0,0].invert_xaxis()

    ax[1,0].plot(rfit,densityfit,color='black')
    ax[1,0].set_ylabel(r'$\rho$')
    ax[1,0].set_yscale('log')
    #ax[1,0].invert_xaxis()

    ax[2,0].plot(rfit,tempfit,color='black')
    ax[2,0].set_ylabel(r'$T$')
    ax[2,0].set_xlabel(r'R [R$_{\odot}$]')

    ax[0,1].plot(rfit,poverrho53fit,color='black')
    #ax[0,1].scatter(q_mass_fraction,poverrho53,marker='x',color='black')
    ax[0,1].set_ylabel(r'${P/{\rho}^(5/3)}$')
    ax[0,1].set_yscale('log')
    #ax[0,1].invert_xaxis()
    
    ax[1,1].plot(rfit,pressurefit,color='black')
    ax[1,1].set_ylabel(r'$P$')
    ax[1,1].set_yscale('log')
    #ax[1,1].invert_xaxis()

    ax[2,1].plot(rfit,jrotfit,color='black')
    ax[2,1].set_ylabel(r'$j$')
    ax[2,1].set_xlabel(r'R [R$_{\odot}$]')
    ax[2,1].set_yscale('log')

    plt.show()

    # graphs H1 and He4

    f,ax=plt.subplots()

    #ax.plot(q,h1fit,color='mediumpurple',label='H1 fraction')
    ax.plot(xq,h1,color='mediumpurple',label='H1 fraction',linewidth=0.75)
    ax.plot(xq,he4,color='blue',label='He4 fraction',linewidth=0.75)
    #ax.plot(q,he4fit,color='blue',label='He4 fraction')
    #ax.plot(q_mass_fraction,h1list,color='black',marker='x')
    #ax.plot(q_mass_fraction,he4list,color='black',marker='x')
    ax.set_title('Star Composition')
    ax.set_xlabel(r'$xq$')
    ax.set_ylabel(r'Fraction')
    ax.invert_xaxis()

    ax.legend(loc='best')
        
    plt.show()


    # graphs H1 and He4 vs radius

    f,ax=plt.subplots()
    
    ax.plot(rcomp,h1,color='mediumpurple',label='H1',linewidth=0.75)
    ax.plot(rcomp,he4,color='blue',label='He4',linewidth=0.75)
    ax.set_title('Star Composition')
    ax.set_xlabel(r'R [R$_{\odot}$]')
    ax.set_ylabel('Fraction')

    plt.show()

    return interp_data
