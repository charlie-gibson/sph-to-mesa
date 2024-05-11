import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
from scipy.optimize import minimize
from eos_func import read_eos
from eos_func import useeostable

def hse_func(starnum=1):

    # Define a function to read EOS data from file or cache
    def read_or_load_eos(eosfile):
        # Check if the cached data exists
        cache_file = "../" + eosfile + ".cache"
        if os.path.exists(cache_file):
            # If the cache file exists, load the data from it
            with open(cache_file, "rb") as f:
                eos_data = pickle.load(f)
        else:
            # Otherwise, read the data from the file and cache it
            eos_data = read_eos(eosfile)
            with open(cache_file, "wb") as f:
                pickle.dump(eos_data, f)
        return eos_data

    uCM = 1
    uGR = 1;
    uSEC = 1
    uERG = 1
    uKELVIN = 1
    pi = np.pi

    uMSUN       = 1.981e33 * uGR
    uRSUN       = 6.9599e10 * uCM
    uG          = 6.67259e-8 * uCM*uCM*uCM/uGR/uSEC**2
    uM_P        = 1.6726231e-24 * uGR
    uC          = 2.99792458e10 * uCM/uSEC
    uK          = 1.380658e-16 * uERG/uKELVIN
    uSIGMA_RAD  = 5.67051e-5 * uERG/uSEC/uCM**2/uKELVIN**4
    uA_RAD      = 4.0*uSIGMA_RAD/uC


    # Initialize empty NumPy arrays
    ammrho = np.array([])
    rsph = np.array([])
    pgassph = np.array([])
    rhosph = np.array([])
    usph = np.array([])
    jrot = np.array([])
    Tsph = np.array([])
    h1 = np.array([])
    he3 = np.array([])
    he4 = np.array([])
    c12 = np.array([])
    n14 = np.array([])
    o16 = np.array([])
    ne20 = np.array([])
    mg24 = np.array([])

    # starnum = int(input("Which star? "))

    # Open the data file for reading
    with open(f'bestfit{starnum}.sph', 'r') as file:
        # Read each line in the file
        for line in file:
            # Split the line into individual values
            values = line.split()
            
            # Convert each value to float and append to respective NumPy arrays
            ammrho = np.append(ammrho, float(values[0]))
            rsph = np.append(rsph, float(values[1]))
            pgassph = np.append(pgassph, float(values[2]))
            rhosph = np.append(rhosph, float(values[3]))
            usph = np.append(usph, float(values[4]))
            jrot = np.append(jrot, float(values[5]))
            Tsph = np.append(Tsph, float(values[6]))
            h1 = np.append(h1, float(values[7]))
            he3 = np.append(he3, float(values[8]))
            he4 = np.append(he4, float(values[9]))
            c12 = np.append(c12, float(values[10]))
            n14 = np.append(n14, float(values[11]))
            o16 = np.append(o16, float(values[12]))
            ne20 = np.append(ne20, float(values[13]))
            mg24 = np.append(mg24, float(values[14]))

    # Print the first few values of A for verification
    #print("beta:", beta[:])
    #print("First few values of A:", A[:])

    def beta_func(beta, y):
        #zeta = 5.0*np.log(1.0 - beta) - 8.0*np.log(beta) + 32.0/beta
        zeta = 5.0*np.log(1.0 - beta) - 8.0*np.log(beta) + 32.0/beta - 32.0
        zeta = y - zeta
        return zeta

    def compute_beta(ptot, A, mean_mu):
        if ptot <= 0:
            return 1

        eps = 1.0e-7
        max_iter = -int(np.log(eps)/np.log(2.0)) + 1

        beta_min = eps
        beta_max = 1.0 - eps

        delta = 3*np.power(uK, 4.0)/uA_RAD

        y = 3.0*np.log(ptot) - 5*np.log(delta) + 12*np.log(A) + 20*np.log(mean_mu*uM_P)

        beta = 0
        for i in range(max_iter):
            beta = 0.5 * (beta_min + beta_max)
            zeta = beta_func(beta, y)
            if zeta < 0:
                beta_min = beta
            else:
                beta_max = beta

        beta = 0.5 * (beta_min + beta_max)
        return beta

    #     H1 He3 He4 C12 N14 O16 Ne20 Mg24
    #Amass = [1, 3, 4, 12, 14, 16, 20, 24]
    #  
    #def am(x):
    #    return (1.0 + Amass[x] / 2.0) / Amass[x]
    #
    #mean_mu = 1 / (
    #        (1.0 + 1)/1 * h1 + \
    #        (1.0 + 2)/3 * he3 + \
    #        am(2) * he4 + \
    #        am(3) * c12 + \
    #        am(4) * n14 + \
    #        am(5) * o16 + \
    #        am(6) * ne20 + \
    #        am(7) * mg24)

    delta = 3*np.power(uK, 4.0)/uA_RAD

    T_from_eos = np.zeros(len(ammrho))
    mean_mu_from_eos = np.zeros(len(ammrho))
    ptot_from_eos = np.zeros(len(ammrho))
    pgas_from_eos = np.zeros(len(ammrho))
    beta_from_eos = np.zeros(len(ammrho))
    A_from_eos = np.zeros(len(ammrho))
    beta_from_solving = np.zeros(len(ammrho))
    rho_from_solving = np.zeros(len(ammrho))

    eos_data = read_or_load_eos('sph.eos_X0.00to0.75step0.05')
    # print(eos_data)

    for index in range(len(ammrho)):
        T_from_eos[index] = useeostable(eos_data=eos_data, ucgs=usph[index], rhocgs=rhosph[index], xxx=h1[index], which=0)
        mean_mu_from_eos[index] = useeostable(eos_data=eos_data, ucgs=usph[index], rhocgs=rhosph[index], xxx=h1[index], which=1)/uM_P
        ptot_from_eos[index] = useeostable(eos_data=eos_data, ucgs=usph[index], rhocgs=rhosph[index], xxx=h1[index], which=2)
        pgas_from_eos[index] = ptot_from_eos[index] - uA_RAD * T_from_eos[index]**4/3
        beta_from_eos[index] = pgas_from_eos[index] / ptot_from_eos[index]

        A_from_eos[index] = pgas_from_eos[index] / rhosph[index] ** (5/3) * np.exp(8/3*(1/beta_from_eos[index] -1))

        beta_from_solving[index] = compute_beta(ptot_from_eos[index], A_from_eos[index], mean_mu_from_eos[index])
        rho_from_solving[index] = delta/A_from_eos[index]**3/(mean_mu_from_eos[index]*uM_P)**4*(1/beta_from_solving[index] -1)*np.exp(8*(1/beta_from_solving[index] - 1))

        #print("beta, beta:", beta[index], result)
        #print("pgas, pgas:", pgassph[index], A[index]*rhosph[index]**(5/3)*np.exp(-8/3*(1/result-1)))
        #print("rho, rho:", rhosph[index], rho2)

    if False:
        # Plots 
        #plt.plot(ammrho/uMSUN, Tsph, 'o', markersize=2,label='T from sph')
        plt.plot(ammrho/uMSUN, T_from_eos, 'o', markersize=2,label='T from EOS')
        plt.xlabel(r'$m\ [M_\odot]$')
        plt.ylabel(r'$T$ [K]')
        plt.grid(True)
        plt.yscale('log')  # Set logarithmic scale on the vertical axis
        #plt.legend()
        plt.show()
        
        plt.plot(ammrho/uMSUN, beta_from_eos, 'o', markersize=2,label=r'from EOS')
        plt.plot(ammrho/uMSUN, beta_from_solving, 'o', markersize=2,label=r'from solving')
        plt.xlabel(r'$m\ [M_\odot]$')
        plt.ylabel(r'$\beta$')
        plt.grid(True)
        #plt.yscale('log')  # Set logarithmic scale on the vertical axis
        plt.legend()
        plt.show()
        
        plt.plot(ammrho/uMSUN, rhosph, 'o', markersize=2,label=r'from bestfit')
        plt.plot(ammrho/uMSUN, rho_from_solving, 'o', markersize=2,label=r'from solving')
        plt.xlabel(r'$m\ [M_\odot]$')
        plt.ylabel(r'$\rho [cgs]$')
        plt.grid(True)
        #plt.yscale('log')  # Set logarithmic scale on the vertical axis
        plt.legend()
        plt.show()
        
        # Create a figure and define subplots
        fig, axs = plt.subplots(5, 1, figsize=(8, 6))
        
        # Plot the first curve in the first subplot
        axs[0].plot(ammrho/uMSUN, T_from_eos, 'o', markersize=2, label='T_from_eos')
        #axs[0].plot(ammrho/uMSUN, Tsph, 'o', markersize=2, label='T from sph')
        axs[0].set_ylabel(r'$T$ [K]')
        axs[0].grid(True)
        axs[0].set_yscale('log')  # Set logarithmic scale on the vertical axis
        
        # Plot the second curve in the second subplot
        axs[1].plot(ammrho/uMSUN, mean_mu_from_eos, 'o', markersize=2, label=r'$\mu$ from EOS')
        #axs[1].plot(ammrho/uMSUN, mean_mu, 'o', markersize=2, label=r'$\mu$')
        axs[1].set_xlabel(r'$m\ [M_\odot]$')
        axs[1].set_ylabel(r'$\mu$')
        axs[1].grid(True)
        axs[1].set_yscale('log')  # Set logarithmic scale on the vertical axis
        
        axs[2].plot(ammrho/uMSUN, pgas_from_eos, 'o', markersize=2, label=r'from EOS')
        axs[2].plot(ammrho/uMSUN, pgassph, 'o', markersize=2, label='from bestfit')
        axs[2].set_xlabel(r'$m\ [M_\odot]$')
        axs[2].set_ylabel(r'$P_{\rm gas}$')
        axs[2].grid(True)
        axs[2].set_yscale('log')  # Set logarithmic scale on the vertical axis
        
        axs[3].plot(ammrho/uMSUN, ptot_from_eos, 'o', markersize=2, label=r'from EOS')
        #axs[3].plot(ammrho/uMSUN, ptot, 'o', markersize=2, label='from bestfit')
        axs[3].set_xlabel(r'$m\ [M_\odot]$')
        axs[3].set_ylabel(r'$P_{\rm tot}$')
        axs[3].grid(True)
        axs[3].set_yscale('log')  # Set logarithmic scale on the vertical axis
        
        axs[4].plot(ammrho/uMSUN, beta_from_eos, 'o', markersize=2, label=r'from EOS')
        #axs[4].plot(ammrho/uMSUN, beta, 'o', markersize=2, label='from bestfit')
        axs[4].set_xlabel(r'$m\ [M_\odot]$')
        axs[4].set_ylabel(r'$\beta$')
        axs[4].grid(True)
        axs[4].set_yscale('log')  # Set logarithmic scale on the vertical axis
        
        # Add legend to each subplot
        axs[0].legend()
        axs[1].legend()
        axs[2].legend()
        axs[3].legend()
        axs[4].legend()
        
        plt.subplots_adjust(hspace=0)
        
        plt.show()


    #Start integration of equation of hydrostatic equilibrium

    Pcmin = 0.9 * ptot_from_eos[0]
    Pcmax = 1e30

    Pctry = Pcmin
    Psurf = ptot_from_eos[-1]

    P = np.zeros(len(ammrho))
    rho = np.zeros(len(ammrho))
    r = np.zeros(len(ammrho))
    beta = np.zeros(len(ammrho))

    converged = False

    # Define the objective function
    def objective_function(u_index):
        # Calculate pressure using useeostable function
        pressure_calculated = useeostable(eos_data=eos_data, ucgs=u_index, rhocgs=rho[index], xxx=h1[index], which=2)
        
        # Calculate the absolute difference between calculated and desired pressure
        difference = abs(pressure_calculated - P[index])
        
        return difference


    while not converged:
        print(f'Trying central pressure {Pctry} [cgs]')
        for index in range(len(ammrho)):
            # Hydrostatic equilibrium requires
            #    dP/dr = - g rho = - G m rho/r^2
            #    dm/dr = 4 pi r^2 rho
            # Because we want to evaluate at the same enclosed mass m values, we'll use m
            # as the integration variable:
            #    dP/dm=dP/dr/(dm/dr)=-G m rho/r^2/(4 pi r^2 rho) = -G m / (4 pi r^4) 
            #    dr/dm = 1/(4 pi r^2 rho)
            
            if index == 0:
                P[0] = Pctry
                # For now, assume mean_mu doesn't change in contraction.  This may need to be fixed.
                beta[0] = compute_beta(P[0], A_from_eos[0], mean_mu_from_eos[0])
                rho[0] = (beta[0]*P[0]/A_from_eos[0] * np.exp(8/3*(1/beta[0] - 1)))**0.6
                # assume constant density out to the the innermost shell: 4/3 pi r^3 = m   =>   r = (3 m/(4 pi))^(1/3)
                r[0] = (3 * ammrho[0]/(4 *pi * rho[0]))**(1/3)
            else:
                dm = ammrho[index] - ammrho[index-1]
                P[index] = P[index-1] - uG/(4*pi) * ammrho[index-1] / r[index-1]**4 * dm
                # For now, assume mean_mu doesn't change in contraction.  This may need to be fixed.
                beta[index] = compute_beta(P[index], A_from_eos[index], mean_mu_from_eos[index])
                rho[index] = (beta[index]*np.abs(P[index])/A_from_eos[index] * np.exp(8/3*(1/beta[index] - 1)))**0.6
                r[index] = r[index-1] + dm/(4*pi*r[index-1]**2*rho[index-1])
                
            # Break the for loop if P[index] becomes negative
            if P[index] <= 0:
                break
        print('Done with for loop at index=',index, ammrho[index], P[index])
        if P[index] <= 0:
            print(f"P[index] became negative when index={index} out of {len(ammrho)-1}")
            Pcmin = Pctry
            if Pcmax >= 1e30:
                Pctry *= 2
            else:
                Pctry = 0.5 * (Pcmin + Pcmax)
        else:
            initial_guess = 1.5 * P[index]/rho[index]
            # Perform optimization to minimize the objective function
            result = minimize(objective_function, initial_guess, tol=1e-8)
            Tsurf = useeostable(eos_data=eos_data, ucgs=result.x[0], rhocgs=rho[index], xxx=h1[index], which=0)
            pgassurf = P[index] - uA_RAD * Tsurf**4/3
            #print('deciding:',P[index],Tsurf,pgassurf,index)
            if pgassurf <= 0:
                Pcmin = Pctry
                if Pcmax >= 1e30:
                    Pctry *= 2
                else:
                    Pctry = 0.5 * (Pcmin + Pcmax)
            elif pgassurf > 2*Psurf:
                Pcmax = Pctry
                Pctry = 0.5 * (Pcmin + Pcmax)
            else:
                converged = True

    u = np.zeros(len(ammrho))
    T = np.zeros(len(ammrho))
    pgas = np.zeros(len(ammrho))
    mean_mu = np.zeros(len(ammrho))

    with open(f"bestfit{starnum}_hse.sph", "w") as f:
        for index in range(len(ammrho)):
            # Perform optimization to minimize the objective function
            initial_guess = 1.5 * P[index]/rho[index]
            result = minimize(objective_function, initial_guess, tol=1e-8)
            # Extract the optimal value of u[index] from the optimization result
            u[index] = result.x[0]        

            # Now u[index] gives the desired pressure
            #print("Optimal value of u[index]:", u[index])
            #print("Compare:", P[index], useeostable(eos_data=eos_data, ucgs=u[index], rhocgs=rho[index], xxx=h1[index], which=2))

            T[index] = useeostable(eos_data=eos_data, ucgs=u[index], rhocgs=rho[index], xxx=h1[index], which=0)
            mean_mu[index] = useeostable(eos_data=eos_data, ucgs=u[index], rhocgs=rho[index], xxx=h1[index], which=1)/uM_P
            pgas[index] = P[index] - uA_RAD * T[index]**4/3

            #if index == 100:
            #    print('about to write:',P[index],T[index],pgas[index])

            f.write('{:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e} {:.9e}\n'.format(
                ammrho[index],
                r[index],
                pgas[index],
                rho[index],
                u[index],
                jrot[index],
                T[index],
                h1[index], he3[index], he4[index], c12[index], n14[index],
                o16[index], ne20[index], mg24[index]))

    if False:
        # Plots 
        
        plt.plot(ammrho/uMSUN, beta_from_eos, 'o', markersize=2,label=r'initial (from eos)')
        plt.plot(ammrho/uMSUN, beta, 'o', markersize=2,label=r'final')
        plt.xlabel(r'$m\ [M_\odot]$')
        plt.ylabel(r'$\beta$')
        plt.grid(True)
        #plt.yscale('log')  # Set logarithmic scale on the vertical axis
        plt.legend()
        plt.show()
        
        plt.plot(ammrho/uMSUN, rhosph, 'o', markersize=2,label=r'sph')
        plt.plot(ammrho/uMSUN, rho, 'o', markersize=2,label=r'final')
        plt.xlabel(r'$m\ [M_\odot]$')
        plt.ylabel(r'$\rho [cgs]$')
        plt.grid(True)
        plt.yscale('log')  # Set logarithmic scale on the vertical axis
        plt.legend()
        plt.show()
        
        plt.plot(ammrho/uMSUN, rsph/uRSUN, 'o', markersize=2,label=r'sph')
        plt.plot(ammrho/uMSUN, r/uRSUN, 'o', markersize=2,label=r'final')
        plt.xlabel(r'$m\ [M_\odot]$')
        plt.ylabel(r'$r\ [R_\odot]$')
        plt.grid(True)
        plt.yscale('log')  # Set logarithmic scale on the vertical axis
        plt.legend()
        plt.show()
        
        # Create a figure and define subplots
        fig, axs = plt.subplots(5, 1, figsize=(8, 6))
        
        # Plot the first curve in the first subplot
        axs[0].plot(ammrho/uMSUN, T_from_eos, 'o', markersize=2, label='original')
        axs[0].plot(ammrho/uMSUN, T, 'o', markersize=2, label='final')
        axs[0].set_ylabel(r'$T$ [K]')
        axs[0].grid(True)
        axs[0].set_yscale('log')  # Set logarithmic scale on the vertical axis
        
        # Plot the second curve in the second subplot
        axs[1].plot(ammrho/uMSUN, mean_mu_from_eos, 'o', markersize=2, label=r'original')
        axs[1].plot(ammrho/uMSUN, mean_mu, 'o', markersize=2, label='final')
        axs[1].set_xlabel(r'$m\ [M_\odot]$')
        axs[1].set_ylabel(r'$\mu$')
        axs[1].grid(True)
        #axs[1].set_yscale('log')  # Set logarithmic scale on the vertical axis
        
        axs[2].plot(ammrho/uMSUN, pgas_from_eos, 'o', markersize=2, label=r'original')
        axs[2].plot(ammrho/uMSUN, pgas, 'o', markersize=2, label='from bestfit')
        axs[2].set_xlabel(r'$m\ [M_\odot]$')
        axs[2].set_ylabel(r'$P_{\rm gas}$')
        axs[2].grid(True)
        axs[2].set_yscale('log')  # Set logarithmic scale on the vertical axis
        
        axs[3].plot(ammrho/uMSUN, ptot_from_eos, 'o', markersize=2, label='original')
        axs[3].plot(ammrho/uMSUN, P, 'o', markersize=2, label='final')
        axs[3].set_xlabel(r'$m\ [M_\odot]$')
        axs[3].set_ylabel(r'$P_{\rm tot}$')
        axs[3].grid(True)
        axs[3].set_yscale('log')  # Set logarithmic scale on the vertical axis
        
        axs[4].plot(ammrho/uMSUN, beta_from_eos, 'o', markersize=2, label='original')
        axs[4].plot(ammrho/uMSUN, beta, 'o', markersize=2, label='final')
        axs[4].set_xlabel(r'$m\ [M_\odot]$')
        axs[4].set_ylabel(r'$\beta$')
        axs[4].grid(True)
        axs[4].set_yscale('log')  # Set logarithmic scale on the vertical axis
        
        # Add legend to each subplot
        axs[0].legend()
        axs[1].legend()
        axs[2].legend()
        axs[3].legend()
        axs[4].legend()
        
        plt.subplots_adjust(hspace=0)
        
        plt.show()


