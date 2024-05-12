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

<<<<<<< HEAD
def entropy_reader(mode,path,starnum=1,opt=0):
=======
<<<<<<< HEAD
def entropy_reader(mode,path,starnum=1,opt=0):
=======
def entropy_reader(mode,path,starnum=1):
>>>>>>> 80b51c09b6d14c8c80b29d98c4f2b29bcd6c776d
>>>>>>> 706dca6185b661b857fdf6d3337257d9d53c74ff

    import numpy as np
    from scipy.interpolate import CubicSpline
    import matplotlib.pyplot as plt
    from eos_func_with_entropy import read_eos
    from eos_func_with_entropy import useeostable

    runit = 6.9599e10

    if opt==0:
        end_string = '.sph'
        end_dat_string = '.dat'
    elif opt==1:
        end_string = '_hse.sph'
        end_dat_string = '_hse.dat'

    # reads in the data from bestfit.sph
<<<<<<< HEAD
    with open(f"bestfit{starnum}"+end_string) as f:
=======
<<<<<<< HEAD
    with open(f"bestfit{starnum}"+end_string) as f:
=======
    with open(f"bestfit{starnum}.sph") as f:
>>>>>>> 80b51c09b6d14c8c80b29d98c4f2b29bcd6c776d
>>>>>>> 706dca6185b661b857fdf6d3337257d9d53c74ff
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
    with open(f"{path}/angular_momentum.dat") as f:
        
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
    
    s = []
    if mode == 'S':
        g = 6.6739e-8 # gravitational constant: cm^3 / g / s^2
        runit = 6.9599e10 # cm
        munit = 1.9891e33 # g
        eosfile = f'{path}/sph.eos_X0.00to0.75step0.05_entropy'
        rhocgs = []
        ucgs = []
        for rho,u in zip(density,specThermEnergy):
            rhocgs.append(rho*munit/runit**3)
            ucgs.append(u*g*munit**2/runit)
        eos_data = read_eos(eosfile)
        for i in range(len(rhocgs)):
            s_val = useeostable(eos_data, ucgs[i], rhocgs[i], h1list[i], which=3)
            s.append(s_val)

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
    if mode=='S':
        s.reverse()

<<<<<<< HEAD
    # uses linear interpolation for all values
=======
<<<<<<< HEAD
    # uses linear interpolation for all values
=======
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

#    xq_mass_fraction.append((xq_mass_fraction[-1]+1)/2)
#    q_mass_fraction.append((q_mass_fraction[-1])/2)
#    mass.append(mass[-1])
#    radius.append(radius[-1])
#    pressure.append(pressure[-1])
#    density.append(density[-1])
#    specThermEnergy.append(specThermEnergy[-1])
#    jrotlist.append(jrotlist[-1])
#    templist.append(templist[-1])
#    h1list.append(h1list[-1])
#    he3list.append(he3list[-1])
#    he4list.append(he4list[-1])
#    c12list.append(c12list[-1])
#    n14list.append(n14list[-1])
#    o16list.append(o16list[-1])
#    ne20list.append(ne20list[-1])
#    mg24list.append(mg24list[-1])
#    poverrho53.append(poverrho53[-1])

    
#    xq_mass_fraction.append(1.0)
#    q_mass_fraction.append(0.0)
#    mass.append(mass[-1])
#    radius.append(radius[-1])
#    pressure.append(pressure[-1])
#    density.append(density[-1])
#    specThermEnergy.append(specThermEnergy[-1])
#    jrotlist.append(jrotlist[-1])
#    templist.append(templist[-1])
#    h1list.append(h1list[-1])
#    he3list.append(he3list[-1])
#    he4list.append(he4list[-1])
#    c12list.append(c12list[-1])
#    n14list.append(n14list[-1])
#    o16list.append(o16list[-1])
#    ne20list.append(ne20list[-1])
#    mg24list.append(mg24list[-1])
#    poverrho53.append(poverrho53[-1])
    
#    xq_mass_fraction.append(1.01)
#    q_mass_fraction.append(-0.01)
#    mass.append(mass[-1])
#    radius.append(radius[-1])
#    pressure.append(pressure[-1])
#    density.append(density[-1])
#    specThermEnergy.append(specThermEnergy[-1])
#    jrotlist.append(jrotlist[-1])
#    templist.append(templist[-1])
#    h1list.append(h1list[-1])
#    he3list.append(he3list[-1])
#    he4list.append(he4list[-1])
#    c12list.append(c12list[-1])
#    n14list.append(n14list[-1])
#    o16list.append(o16list[-1])
#    ne20list.append(ne20list[-1])
#    mg24list.append(mg24list[-1])
#    poverrho53.append(poverrho53[-1])

#    xq_mass_fraction.append(1.02)
#    q_mass_fraction.append(-0.02)
#    mass.append(mass[-1])
#    radius.append(radius[-1])
#    pressure.append(pressure[-1])
#    density.append(density[-1])
#    specThermEnergy.append(specThermEnergy[-1])
#    jrotlist.append(jrotlist[-1])
#    templist.append(templist[-1])
#    h1list.append(h1list[-1])
#    he3list.append(he3list[-1])
#    he4list.append(he4list[-1])
#    c12list.append(c12list[-1])
#    n14list.append(n14list[-1])
#    o16list.append(o16list[-1])
#    ne20list.append(ne20list[-1])
#    mg24list.append(mg24list[-1])
#    poverrho53.append(poverrho53[-1])

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
>>>>>>> 80b51c09b6d14c8c80b29d98c4f2b29bcd6c776d
>>>>>>> 706dca6185b661b857fdf6d3337257d9d53c74ff
    h1 = np.interp(x=xq,xp=xq_mass_fraction, fp=h1list)
    he3 = np.interp(x=xq,xp=xq_mass_fraction, fp=he3list)
    he4 = np.interp(x=xq,xp=xq_mass_fraction, fp=he4list)
    c12 = np.interp(x=xq,xp=xq_mass_fraction, fp=c12list)
    n14 = np.interp(x=xq,xp=xq_mass_fraction, fp=n14list)
    o16 = np.interp(x=xq,xp=xq_mass_fraction, fp=o16list)
    ne20 = np.interp(x=xq,xp=xq_mass_fraction, fp=ne20list)
    mg24 = np.interp(x=xq,xp=xq_mass_fraction, fp=mg24list)
    rcomp = np.interp(x=xq,xp=xq_mass_fraction, fp=radius)
    jrot = np.interp(x=xq,xp=xq_mass_fraction, fp=jrotlist)
    rhofit = np.interp(x=xq,xp=xq_mass_fraction, fp=density)
    pressurefit = np.interp(x=xq,xp=xq_mass_fraction, fp=pressure)
    specThermEnergyfit = np.interp(x=xq,xp=xq_mass_fraction, fp=specThermEnergy)
    poverrho53fit = np.interp(x=xq,xp=xq_mass_fraction, fp=poverrho53)
    tempfit = np.interp(x=xq,xp=xq_mass_fraction, fp=templist)
    rfit = np.interp(x=xq,xp=xq_mass_fraction, fp=radius)
    if mode=='S':
        sfit = np.interp(x=xq,xp=xq_mass_fraction,fp=s)

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
    if mode=='S':
        s.reverse()


    # creates the best fit equations to be graphed later
    # densityfit = rho(xq)
    # rfit = r(xq)
    # specThermEnergyfit = e(xq)
    # poverrho53fit = prho53(xq)
    # pressurefit = p(xq)
    # jrotfit = jrot(xq)
    # tempfit = temp(xq)
    #h1fit = h1(xq)
    #he3fit = he3(xq)
    #he4fit = he4(xq)
    #c12fit = c12(xq)
    #n14fit = n14(xq)
    #o16fit = o16(xq)
    #ne20fit = ne20(xq)
    #mg24fit = mg24(xq)

    try:
        interp_data = {
            'xq':xq,
            'q':q,
            'r':rfit,
            'rho':rhofit,
            'P':pressurefit,
            'A':poverrho53fit,
            'T':tempfit,
            'jrot':jrot,
            'H1':h1,
            'He3':he3,
            'He4':he4,
            'C12':c12,
            'N14':n14,
            'O16':o16,
            'Ne20':ne20,
            'Mg24':mg24,
            'S':s
        }
    except:        
        interp_data = {
            'xq':xq,
            'q':q,
            'r':rfit,
            'rho':rhofit,
            'P':pressurefit,
            'A':poverrho53fit,
            'T':tempfit,
            'jrot':jrot,
            'H1':h1,
            'He3':he3,
            'He4':he4,
            'C12':c12,
            'N14':n14,
            'O16':o16,
            'Ne20':ne20,
            'Mg24':mg24
        }
    
<<<<<<< HEAD
    mycompositionfile = open(f"composition{starnum}"+end_dat_string, "w")
=======
<<<<<<< HEAD
    mycompositionfile = open(f"composition{starnum}"+end_dat_string, "w")
=======
    mycompositionfile = open(f"composition{starnum}.dat", "w")
>>>>>>> 80b51c09b6d14c8c80b29d98c4f2b29bcd6c776d
>>>>>>> 706dca6185b661b857fdf6d3337257d9d53c74ff

    mycompositionfile.write(f"{num_zones[0]}    {8}\n")
        
    for k in range(len(xq)):
        mycompositionfile.write(f"{xq[k]}     {h1[k]}     {he3[k]}     {he4[k]}     {c12[k]}     {n14[k]}     {o16[k]}     {ne20[k]}     {mg24[k]}\n")
        

    # this writes the entropy.dat file to be used in mesa

<<<<<<< HEAD
    myentropyfile = open(f"entropy{starnum}"+end_dat_string, "w")
=======
<<<<<<< HEAD
    myentropyfile = open(f"entropy{starnum}"+end_dat_string, "w")
=======
    myentropyfile = open(f"entropy{starnum}.dat", "w")
>>>>>>> 80b51c09b6d14c8c80b29d98c4f2b29bcd6c776d
>>>>>>> 706dca6185b661b857fdf6d3337257d9d53c74ff

    # includes the number of zones at the top of the file
    myentropyfile.write(f"{num_zones[0]}\n")

    if mode == 'DT':

        # this format is for a density, temperature pair
        for j in range(len(xq)):
            myentropyfile.write(f"{xq[j]}     {rhofit[j]}     {tempfit[j]}\n")

        myentropyfile.close()
    
    elif mode == 'PT':

        # print('length of q =', len(q))
        for j in range(len(xq)):
            myentropyfile.write(f"{xq[j]}     {pressurefit[j]}     {tempfit[j]}\n")
        
        myentropyfile.close

    elif mode == 'DE':

        # this format is for a density, specific thermal energy pair
        # writes mass fraction, density, and specific thermal energy on one line per shell
        for j in range(len(xq)):
            myentropyfile.write(f"{xq[j]}     {rhofit[j]}     {specThermEnergyfit[j]}\n")

        myentropyfile.close
    
    elif mode == 'S':
        for j in range(len(xq)):
            myentropyfile.write(f'{xq[j]}     {sfit[j]}\n')

        myentropyfile.close()

    elif mode == 'DP':

        for j in range(len(xq)):
            myentropyfile.write(f"{j}     {rhofit[j]}     {pressurefit[j]}\n")
        
        myentropyfile.close()

    # this writes the angular_momentum.dat file to be used in mesa

<<<<<<< HEAD
    myjrotfile = open(f"angular_momentum{starnum}"+end_dat_string, "w")
=======
<<<<<<< HEAD
    myjrotfile = open(f"angular_momentum{starnum}"+end_dat_string, "w")
=======
    myjrotfile = open(f"angular_momentum{starnum}.dat", "w")
>>>>>>> 80b51c09b6d14c8c80b29d98c4f2b29bcd6c776d
>>>>>>> 706dca6185b661b857fdf6d3337257d9d53c74ff

    # includes the number of zones at the top of the file
    myjrotfile.write(f"{num_zones[0]}\n")

    # writes mass fraction, density, and specific thermal energy on one line per shell
    for k in range(len(xq)):
        myjrotfile.write(f"{xq[k]}     {jrot[k]}\n")

    myjrotfile.close()

    # calculates the total X, Y, and Z

    weightedH1=0
    weightedHe3=0
    weightedHe4=0
    weightedC12=0
    weightedN14=0
    weightedO16=0
    weightedNe20=0
    weightedMg24=0

    for i in range(len(xq)-1):
        xq_frac=xq[i+1]-xq[i]
        weightedH1+=h1[i]*xq_frac
        weightedHe3+=he3[i]*xq_frac
        weightedHe4+=he4[i]*xq_frac
        weightedC12+=c12[i]*xq_frac
        weightedN14+=n14[i]*xq_frac
        weightedO16+=o16[i]*xq_frac
        weightedNe20+=ne20[i]*xq_frac
        weightedMg24+=mg24[i]*xq_frac

    Z=1-weightedH1-weightedHe4
    print(f'X: {weightedH1}\nY: {weightedHe4}\nZ: {Z}\n')
    print(f'H1: {weightedH1}\nHe3: {weightedHe3}\nHe4: {weightedHe4}\nC12: {weightedC12}\nN14: {weightedN14} \nO16: {weightedO16}\nNe20: {weightedNe20}\nMg24: {weightedMg24}')

<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> 706dca6185b661b857fdf6d3337257d9d53c74ff
    with open(f'sph_star{starnum}.dat','r') as f:
        # Mass:                              0.890972185141951
        # Radius:                            89.41858590668664
        # Particles:                         10421
        # Particles for inner 50% of mass:   3054
        # Star 1 Particles:                  10421
        # Star 2 Particles:                  0
        # Star 1 Fractional Mass:            1.0000000000000067
        # Star 2 Fractional Mass:            0.0
        # Total Angular Momentum:            6.313800040147736e+49
        # Total Specific Angular Momentum:   3.562625e+16
        n=0
        for line in f:
            val=float(line.split(':')[1])
            # print(val)
            if n==0:
                Mass = val
            if n==1:
                if opt == 1:
                    Radius = radius[-1] / runit
                else:
                    Radius = val
            if n==2:
                Particles = val
            if n==3:
                Particles_50 = val
            if n==4:
                Particles1 = val
            if n==5:
                Particles2 = val
            if n==6:
                Mass1 = val
            if n==7:
                Mass2 = val
            if n==8:
                L = val
            if n==9:
                J = val
            if n>=10:
                break
            n += 1

    with open(f'sph_star{starnum}'+end_dat_string,'w') as f:
            f.write(f'Mass:                             {Mass}\n')
            f.write(f'Radius:                           {Radius}\n')
            f.write(f'Particles:                        {int(Particles)}\n')
            f.write(f'Particles for inner 50% of mass:  {int(Particles_50)}\n')
            f.write(f'Star 1 Particles:                 {int(Particles1)}\n')
            f.write(f'Star 2 Particles:                 {int(Particles2)}\n')
            f.write(f'Fraction of Mass From Star 1:     {Mass1}\n')
            f.write(f'Fraction of Mass From Star 2:     {Mass2}\n')
            f.write(f'Total Angular Momentum:           {L}\n')
            f.write(f'Total Specific Angular Momentum:  {J}\n')
            f.write(f'X:                                {weightedH1}\n')
            f.write(f'Y:                                {weightedHe4}\n')
            f.write(f'Z:                                {Z}\n')
            f.write(f'H1:                               {weightedH1}\n')
            f.write(f'He3:                              {weightedHe3}\n')
            f.write(f'He4:                              {weightedHe4}\n')
            f.write(f'C12:                              {weightedC12}\n')
            f.write(f'N14:                              {weightedN14}\n')
            f.write(f'O16:                              {weightedO16}\n')
            f.write(f'Ne20:                             {weightedNe20}\n')
            f.write(f'Mg24:                             {weightedMg24}\n')
            f.write(f'[N/C]:                            {weightedN14/weightedC12}\n')
<<<<<<< HEAD
=======
=======
    with open(f'sph_star{starnum}.dat','a') as f:
        f.write(f'X:                                {weightedH1}\n')
        f.write(f'Y:                                {weightedHe4}\n')
        f.write(f'Z:                                {Z}\n')
        f.write(f'H1:                               {weightedH1}\n')
        f.write(f'He3:                              {weightedHe3}\n')
        f.write(f'He4:                              {weightedHe4}\n')
        f.write(f'C12:                              {weightedC12}\n')
        f.write(f'N14:                              {weightedN14}\n')
        f.write(f'O16:                              {weightedO16}\n')
        f.write(f'Ne20:                             {weightedNe20}\n')
        f.write(f'Mg24:                             {weightedMg24}\n')
        f.write(f'[N/C]:                            {weightedN14/weightedC12}')
    
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

    #plt.show()

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

    #plt.show()

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
        
    #plt.show()


    # graphs H1 and He4 vs radius

    f,ax=plt.subplots()
    
    ax.plot(rcomp,h1,color='mediumpurple',label='H1',linewidth=0.75)
    ax.plot(rcomp,he4,color='blue',label='He4',linewidth=0.75)
    ax.set_title('Star Composition')
    ax.set_xlabel(r'R [R$_{\odot}$]')
    ax.set_ylabel('Fraction')

    #plt.show()
>>>>>>> 80b51c09b6d14c8c80b29d98c4f2b29bcd6c776d
>>>>>>> 706dca6185b661b857fdf6d3337257d9d53c74ff

    return interp_data
