"""
This code calculates the eccentricity of two stars in an orbit.
It finds the center of each star (based on the densest particle)
and uses the relative velocity and distance between the two points
in the eccentricity equation. It also finds the relative velocity
at infinity. These can be used to find the semi-major axis and other
quantities.

Charles Gibson
Allegheny College
Department of Physics
03/26/2024
"""

import numpy as np
import matplotlib.pyplot as plt
from component_reader import component_reader

# constants
G = 6.6739e-8 # gravitational constant: cm^3 / g / s^2
runit = 6.9599e10 # cm
munit = 1.9891e33 # g
tunit = np.sqrt(runit**3/(G*munit))

# function for eccentricity calcs
# needs the data from the output file and the number of the file
def eccentricity(readit_data, component_file_number):
    # collects the important data from readit
    mass = readit_data['am']
    x = readit_data['x']
    y = readit_data['y']
    z = readit_data['z']
    vx = readit_data['vx']
    vy = readit_data['vy']
    vz = readit_data['vz']
    rho = readit_data['rho']

    # finds the components from the comp*.sph file
    components = component_reader(component_file_number)
    print(f'READ IN comp{component_file_number}.sph')

    # insantiates variables
    m1 = 0
    m2 = 0
    rho1 = 0
    rho2 = 0

    # loops through every particle
    for i in range(len(components)):
        
        # if the particle belongs to star 1, the mass is added to m1
        if components[i]==1:
            m1 += mass[i]
            if rho[i] > rho1:
                # if the particle is the densest particle in star 1
                # it is considered to be the center so the position
                # and velocity of the center is found
                rho1 = rho[i]
                x1 = x[i]
                y1 = y[i]
                z1 = z[i]
                vx1 = vx[i]
                vy1 = vy[i]
                vz1 = vz[i]
        # does the same kind of calcs for star 2
        elif components[i]==2:
            m2 += mass[i]
            if rho[i] > rho2:
                rho2 = rho[i]
                x2 = x[i]
                y2 = y[i]
                z2 = z[i]
                vx2 = vx[i]
                vy2 = vy[i]
                vz2 = vz[i]
    
    # converts positions, velocities, and masses to cgs
    x1 = np.array(x1) * runit
    y1 = np.array(y1) * runit
    z1 = np.array(z1) * runit
    x2 = np.array(x2) * runit
    y2 = np.array(y2) * runit
    z2 = np.array(z2) * runit

    vx1 = np.array(vx1) * runit / tunit
    vy1 = np.array(vy1) * runit / tunit
    vz1 = np.array(vz1) * runit / tunit

    vx2 = np.array(vx2) * runit / tunit
    vy2 = np.array(vy2) * runit / tunit
    vz2 = np.array(vz2) * runit / tunit

    m1 = m1 * munit
    m2 = m2 * munit
    print('m1 and m2: ', m1/munit, m2/munit)

    # calculates the relative velocity and position between the center
    # of the two stars
    r_rel = np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
    v_rel = np.sqrt((vx1-vx2)**2+(vy1-vy2)**2+(vz1-vz2)**2)
    print('Relative r: ',r_rel/runit)
    print('Relative v: ',v_rel/1e5)

    # calculates the important quantities for the eccentricity equation
    mu = m1*m2/(m1+m2)
    k = G*m1*m2
    # l2 = mu**2*r_rel**2*v_rel**2
    l2 = mu**2*np.linalg.norm(np.cross([(x1-x2), (y1-y2), (z1-z2)],[(vx1-vx2), (vy1-vy2), (vz1-vz2)]))**2
    k_energy = 1/2*mu*v_rel**2
    u_energy = -k/r_rel
    # energy = 1/2*mu*v_rel**2-k/r_rel
    energy = k_energy + u_energy

    print('Kinetic Energy: ', k_energy)
    print('Potential Energy: ',u_energy)
    print('Energy: ',energy)
    print('Angular Momentum: ', np.sqrt(l2))
    print('Angular Momentum squared: ',l2)
    print('Mu: ',mu)
    print('k: ',k)

    # calculates eccentricity and relative velocity at infinity
    e = np.sqrt(1+2*energy*l2/(mu*k**2))
    v_inf_rel = np.sqrt(2*energy/mu) / 1e5
    v_inf_rel2 = 2*energy/mu / 1e10
    print(f'Eccentricity: {e}')
    print(f'v_rel_infty: {v_inf_rel}')

    # writes eccentricity and relative velocity at infinity to eccentricity.out
    file = 'eccentricity.out'
    with open(file,'w') as f:
        f.write(f'{e}\n')
#        f.write(f'{v_inf_rel2}\n')
        f.write(f'{v_inf_rel}\n')
