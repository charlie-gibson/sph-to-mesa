def v():
    import matplotlib.pyplot as plt
    import numpy as np


    # reads in energy0.sph: change this if you have higher numbers of simulations in the same directory
    # e.g. energy1.sph, energy2.sph, etc.
    filename = 'energy0.sph'

    # collects the data outputted to energy*.sph
    time = []
    w = []
    t = []
    e = []
    u = []
    s = []
    am = []

    with open(filename, 'r') as file:
        for line in file:
            line=line.split()
            time.append(float(line[0]))
            w.append(float(line[1]))
            t.append(float(line[2]))
            u.append(float(line[3]))
            e.append(float(line[4]))
            s.append(float(line[5]))
            am.append(float(line[6]))

    # creates a set of subplots similar to the original v/vln plots in supermongo
    f,[ax5,ax4,ax3,ax2]=plt.subplots(4,1,sharex=True)
    f.subplots_adjust(hspace=0)

    ax5.plot(time,u,c='black',linewidth=0.75)
    ax5.set_ylabel(r'${U}$')

    ax4.plot(time,w,c='black',linewidth=0.75)
    ax4.set_ylabel(r'${W}$')

    ax3.plot(time,t,c='black',linewidth=0.75)
    ax3.set_ylabel(r'${T}$')

    ax2.plot(time,e,c='black',linewidth=0.75)
    ax2.set_ylabel(r'${E}$')

    # creates another set of subplots for angular momentum and entropy vs time
    f,[ax1,ax0]=plt.subplots(2,1,sharex=True)
    f.subplots_adjust(hspace=0)

    ax1.plot(time,am,c='black',linewidth=0.25)
    ax1.set_ylabel(r'${j}$')

    ax0.plot(time,s,c='black',linewidth=0.75)
    ax0.set_xlabel(r'${t}$')
    ax0.set_ylabel(r'${s}$')

    # outputs both plots at the same time
    plt.show()
