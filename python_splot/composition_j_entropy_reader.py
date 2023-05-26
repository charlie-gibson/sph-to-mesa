"""
This code analyzes the file bestfit.sph to generate
a profile pair of density and specific thermal energy.

It now also reads in data for the angular momentum and composition.

It uses these data to write three files:
composition.dat, entropy.dat, angular_momentum.dat

These values can be used to generate a mesa run using
the relax_composition_j_entropy test suite.

Charles Gibson
Allegheny College
Department of Physics
"""

def entropy_reader():

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
                h1list.append(float(values[6].replace("D", "e")))
                he3list.append(float(values[7].replace("D", "e")))
                he4list.append(float(values[8].replace("D", "e")))
                c12list.append(float(values[9].replace("D", "e")))
                n14list.append(float(values[10].replace("D", "e")))
                o16list.append(float(values[11].replace("D", "e")))
                ne20list.append(float(values[12].replace("D", "e")))
                mg24list.append(float(values[13].replace("D", "e")))
            except:
                pass

    # temporary until we can figure out how num_zones is calculated in MESA
    with open("/home/gibson01/splot_directories/python_splot/angular_momentum.dat") as f:
        
        q = []

        n = 0

        for line in f:
            if n == 0:
                num_zones = line.split()
                n += 1
            else:
                newq = line.split()[0]
                q.append(float(newq))

    myentropyfile = open("entropy.dat", "w")
    totalmass = 0
    q_mass_fraction = []

    # the total mass should be the final line's mass
    # which we find with this line

    totalmass = float(mass[-1])

    index = 0
    # calculates the total mass outside each "shell"
    for n in mass:
        massenclosed = n # each value in the mass list is how much mass is found from the center to that radius
        massoutside = totalmass - massenclosed
        massfraction = massoutside / totalmass
        q_mass_fraction.append(massfraction) # mesa uses this version of the mass fraction
        # print(q_mass_fraction[index])
        index += 1

    # calculates P/rho^(5/3)
    poverrho53 = []
    for i, j in zip(density, pressure):
        poverrho53.append(i/j**(5/3))

    # Spline can only interpolate with increasing values of x
    q_mass_fraction.reverse()
    density.reverse()
    specThermEnergy.reverse()
    poverrho53.reverse()
    jrotlist.reverse()
    # h1list.reverse()
    # he3list.reverse()
    # he4list.reverse()
    # c12list.reverse()
    # n14list.reverse()
    # o16list.reverse()
    # ne20list.reverse()
    # mg24list.reverse()

    # creates an equation for the data of mass fraction vs density,
    # mass fraction vs specific thermal energy, and mass fraction and angular momentum
    rho = CubicSpline(q_mass_fraction, density)
    e = CubicSpline(q_mass_fraction, specThermEnergy)
    prho53 = CubicSpline(q_mass_fraction, poverrho53)
    jrot = CubicSpline(q_mass_fraction, jrotlist)
    h1 = CubicSpline(q_mass_fraction, h1list)
    he3 = CubicSpline(q_mass_fraction, he3list)
    he4 = CubicSpline(q_mass_fraction, he4list)
    c12 = CubicSpline(q_mass_fraction, c12list)
    n14 = CubicSpline(q_mass_fraction, n14list)
    o16 = CubicSpline(q_mass_fraction, o16list)
    ne20 = CubicSpline(q_mass_fraction, ne20list)
    mg24 = CubicSpline(q_mass_fraction, mg24list)

    # now that spline is made, we need to write to the file
    q_mass_fraction.reverse()
    density.reverse()
    specThermEnergy.reverse()
    poverrho53.reverse()
    jrotlist.reverse()
    h1list.reverse()
    he3list.reverse()
    he4list.reverse()
    c12list.reverse()
    n14list.reverse()
    o16list.reverse()
    ne20list.reverse()
    mg24list.reverse()

    # creates the best fit equations to be graphed later
    densityfit = rho(q)
    specThermEnergyfit = e(q)
    poverrho53fit = prho53(q)
    jrotfit = jrot(q)
    h1fit = h1(q)
    he3fit = he3(q)
    he4fit = he4(q)
    c12fit = c12(q)
    n14fit = n14(q)
    o16fit = o16(q)
    ne20fit = ne20(q)
    mg24fit = mg24(q)


    mycompositionfile = open("composition.dat", "w")

    mycompositionfile.write(f"{num_zones[len(h1)]}    {8}\n")

    for k in q:
        mycompositionfile.write(f"{k}     {h1(k)}     {he3(k)}     {he4(k)}     {c12(k)}     \
{n14(k)}     {o16(k)}     {ne20(k)}     {mg24(k)}\n")

    # this writes the entropy.dat file to be used in mesa

    # includes the number of zones at the top of the file
    myentropyfile.write(f"{num_zones[0]}\n")

    # writes mass fraction, density, and specific thermal energy on one line per shell
    for j in q:
        myentropyfile.write(f"{j}     {rho(j)}     {e(j)}\n")

    myentropyfile.close()

    # this writes the angular_momentum.dat file to be used in mesa

    myjrotfile = open("angular_momentum.dat", "w")

    # includes the number of zones at the top of the file
    myjrotfile.write(f"{num_zones[0]}\n")

    # writes mass fraction, density, and specific thermal energy on one line per shell
    for j in q:
        myjrotfile.write(f"{j}     {jrot(j)}\n")

    myjrotfile.close()

    # graphs the data for visual representation
    plt.subplot(2, 2, 1)
    plt.plot(q, specThermEnergyfit)
    plt.title("Specific Thermal Energy", fontsize=9)
    plt.xlabel("$q$")
    plt.ylabel("$e$")
    plt.gca().invert_xaxis()

    plt.subplot(2, 2, 2)
    plt.title("Density", fontsize=9)
    plt.plot(q, densityfit)
    plt.xlabel("$q$")
    plt.ylabel(r"$\rho$")
    plt.gca().invert_xaxis()

    plt.subplot(2, 2, 3)
    plt.title(r"$\frac{P}{\rho^\dfrac{5}{3}}$", fontsize=9)
    plt.plot(q, poverrho53fit)
    plt.xlabel("q")
    plt.ylabel("$A$")
    plt.gca().invert_xaxis()

    plt.subplot(2, 2, 4)
    plt.title("Angular momentum", fontsize=9)
    plt.plot(q, jrotfit)
    plt.xlabel("$q$")
    plt.ylabel("$Angular Momentum$")
    plt.gca().invert_xaxis()

    plt.subplots_adjust(wspace=0.4, hspace=0.4)

    plt.show()
