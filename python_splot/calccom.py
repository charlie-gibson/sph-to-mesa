"""
James Lombardi
"""

def calccom(ntot,am,x,y,z,vx,vy,vz,icomp):
    import numpy as np
    am1 = 0.0
    x1 = 0.0
    y1 = 0.0
    z1 = 0.0
    vx1 = 0.0
    vy1 = 0.0
    vz1 = 0.0
    am2 = 0.0
    x2 = 0.0
    y2 = 0.0
    z2 = 0.0
    vx2 = 0.0
    vy2 = 0.0
    vz2 = 0.0
    am3 = 0.0
    x3 = 0.0
    y3 = 0.0
    z3 = 0.0
    vx3 = 0.0
    vy3 = 0.0
    vz3 = 0.0
    am4 = 0.0

    for i in range(ntot):
        if icomp[i] == 1:
            am1 += am[i]
            x1 += am[i] * x[i]
            y1 += am[i] * y[i]
            z1 += am[i] * z[i]
            vx1 += am[i] * vx[i]
            vy1 += am[i] * vy[i]
            vz1 += am[i] * vz[i]
        elif icomp[i] == 2:
            am2 += am[i]
            x2 += am[i] * x[i]
            y2 += am[i] * y[i]
            z2 += am[i] * z[i]
            vx2 += am[i] * vx[i]
            vy2 += am[i] * vy[i]
            vz2 += am[i] * vz[i]
        elif icomp[i] == 3:
            am3 += am[i]
            x3 += am[i] * x[i]
            y3 += am[i] * y[i]
            z3 += am[i] * z[i]
            vx3 += am[i] * vx[i]
            vy3 += am[i] * vy[i]
            vz3 += am[i] * vz[i]
        elif icomp[i] == 4:
            am4 += am[i]

    if am1 > 0.0:
        x1 /= am1
        y1 /= am1
        z1 /= am1
        vx1 /= am1
        vy1 /= am1
        vz1 /= am1

    if am2 > 0.0:
        x2 /= am2
        y2 /= am2
        z2 /= am2
        vx2 /= am2
        vy2 /= am2
        vz2 /= am2

    if am3 > 0.0:
        x3 /= am3
        y3 /= am3
        z3 /= am3
        vx3 /= am3
        vy3 /= am3
        vz3 /= am3

    return am1, np.array([x1, y1, z1]), np.array([vx1, vy1, vz1]),\
        am2, np.array([x2, y2, z2]), np.array([vx2, vy2, vz2]),\
        am3, np.array([x3, y3, z3]), np.array([vx3, vy3, vz3]),\
        am4
