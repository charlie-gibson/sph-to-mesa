"""
Writes out????.sph to an ascii file
out????_ascii.sph

Charles Gibson
Allegheny College
Department of Physics
05/04/2024
"""

import numpy as np

def ascii_output(nnit, readit_data, header_data):

    if nnit < 10:
        FNAME = f'out000{nnit}_ascii.sph'
    elif nnit < 100:
        FNAME = f'out00{nnit}_ascii.sph'
    elif nnit < 1000:
        FNAME = f'out0{nnit}_ascii.sph'
    else:
        FNAME = f'out{nnit}_ascii.sph'

    print(f'WRITING TO {FNAME}')
    
    ntot = header_data[0]
    nnopt = header_data[1]
    hco = header_data[2]
    hfloor = header_data[3]
    sep0 = header_data[4]
    tf = header_data[5]
    dtout = header_data[6]
    nout = header_data[7]
    nit = header_data[8]
    t = header_data[9]
    nav = header_data[10]
    alpha = header_data[11]
    beta = header_data[12]
    tjumpahead = header_data[13]
    ngr = header_data[14]
    nrelax = header_data[15]
    trelax = header_data[16]
    dt = header_data[17]
    omega2 = header_data[18]
    ncooling = header_data[19]
    erad = header_data[20]
    ndisplace = header_data[21]
    displacex = header_data[22]
    displacey = header_data[23]
    displacez = header_data[24]

    x = readit_data['x']
    y = readit_data['y']
    z = readit_data['z']
    am = readit_data['am']
    hp = readit_data['hp']
    rho = readit_data['rho']
    vx = readit_data['vx']
    vy = readit_data['vy']
    vz = readit_data['vz']
    vxdot = readit_data['vxdot']
    vydot = readit_data['vydot']
    vzdot = readit_data['vzdot']
    u = readit_data['u']
    udot = readit_data['udot']
    grpot = readit_data['grpot']
    meanmolecular = readit_data['meanmolecular']
    cc = readit_data['cc']
    divv = readit_data['divv']

    with open(FNAME,'w') as f:
        f.write('      ntot     nnopt                      hco                   hfloor                     sep0\
                       tf                    dtout      nout       nit                        t       nav                    alpha\
                     beta               tjumpahead       ngr    nrelax                   trelax                       dt\
                   omega2  ncooling                     erad      ndisplace                displacex                displacey\
                displacez\n')
        f.write(f'{ntot:10}{nnopt:10}{hco:25}{hfloor:25}{sep0:25}{tf:25}{dtout:25}{nout:10}{nit:10}{t:25}{nav:10}{alpha:25}{beta:25}{tjumpahead:25}{ngr:10}\
{nrelax:10}{trelax:25}{dt:25}{omega2:25}{ncooling:10}{erad:25}{ndisplace:15}{displacex:25}{displacey:25}{displacez:25}\n')
        f.write('                        x                        y                        z                       am\
                       hp                      rho                       vx                       vy                       vz\
                    vxdot                    vydot                    vzdot                        u                     udot\
                    grpot     meanmolecular[grams]        cc                     divv\n')
        for i in range(ntot):
            f.write(f'{x[i]:25}{y[i]:25}{ z[i]:25}{am[i]:25}{hp[i]:25}{rho[i]:25}{vx[i]:25}{vy[i]:25}{vz[i]:25}{vxdot[i]:25}\
{vydot[i]:25}{vzdot[i]:25}{u[i]:25}{udot[i]:25}{grpot[i]:25}{meanmolecular[i]:25}{cc[i]:10}{divv[i]:25}\n')
        f.write(f'{ntot:10}\n')
