"""
This code uses two files to determine position and component values for
particles that are still bound to a star after an interaction

Charles Gibson
Allegheny College
Department of Physics
"""

def bound_particle_data(data, component, composition,jrot=0,temp=0,press=0,component_val=1):

    import numpy as np

    # uses the lists in data (a dictionary) from readit.py to determine the variables for calculations
    ntot = int(data['ntot']) # number of particles
    x = data['x'] # x position
    y = data['y'] # y position
    z = data['z'] # z position
    am = data['am'] # mass
    hp = data['hp'] # smoothing length
    rho = data['rho'] # density
    vx = data['vx'] # x velocity
    vy = data['vy'] # y velocity
    vz = data['vz'] # z velocity
    vxdot = data['vxdot'] # x acceleration
    vydot = data['vydot'] # y acceleration
    vzdot = data['vzdot'] # z acceleration
    u = data['u'] # specific internal energy
    udot = data['udot'] # d/dt[u]
    grpot = data['grpot'] # gravitational potential
    meanmolecular = data['meanmolecular'] # mean molecular weight
    cc = data['cc']
    divv = data['divv']
    # aa = data['aa']
    # bb = data['bb']
    # dd = data['dd']
    try:
        n1 = data['n1']
        n2 = data['n2']
        cc1val = data['cc1val']
        cc2val = data['cc2val']
    except:
        pass

    h1 = composition[0]
    he3 = composition[1]
    he4 = composition[2]
    c12 = composition[3]
    n14 = composition[4]
    o16 = composition[5]
    ne20 = composition[6]
    mg24 = composition[7]

    #try:
    #    junk=jrot[ntot]
    #except:
    #    jrot=np.zeros(ntot)
    #    press=np.zeros(ntot)
    #    temp=np.zeros(ntot)

    boundx = []
    boundy = []
    boundz = []
    boundam = []
    boundhp = []
    boundrho = []
    boundvx = []
    boundvy = []
    boundvz = []
    boundvxdot = []
    boundvydot = []
    boundvzdot = []
    boundu = []
    boundudot = []
    boundgrpot = []
    boundmeanmolecular = []
    boundcc = []
    bounddivv = []
    bound_jrot = []
    bound_temp = []
    bound_pressure = []
    # boundaa = []
    # boundbb = []
    # bounddd = []
    boundh1 = []
    boundhe3 = []
    boundhe4 = []
    boundc12 = []
    boundn14 = []
    boundo16 = []
    boundne20 = []
    boundmg24 = []

    n = 0

    for value in component:
        if value == component_val:
            boundx.append(x[n])
            boundy.append(y[n])
            boundz.append(z[n])
            boundam.append(am[n])
            boundhp.append(hp[n])
            boundrho.append(rho[n])
            boundvx.append(vx[n])
            boundvy.append(vy[n])
            boundvz.append(vz[n])
            boundvxdot.append(vxdot[n])
            boundvydot.append(vydot[n])
            boundvzdot.append(vzdot[n])
            boundu.append(u[n])
            boundudot.append(udot[n])
            boundgrpot.append(grpot[n])
            boundmeanmolecular.append(meanmolecular[n])
            boundcc.append(cc[n])
            bounddivv.append(divv[n])
            bound_jrot.append(jrot[n])
            bound_temp.append(temp[n])
            bound_pressure.append(press[n])
            # boundaa.append(aa[n])
            # boundbb.append(bb[n])
            # bounddd.append(dd[n])

            boundh1.append(h1[n])
            boundhe3.append(he3[n])
            boundhe4.append(he4[n])
            boundc12.append(c12[n])
            boundn14.append(n14[n])
            boundo16.append(o16[n])
            boundne20.append(ne20[n])
            boundmg24.append(mg24[n])
        
        n += 1
    
    boundx = np.array(boundx)
    boundy = np.array(boundy)
    boundz = np.array(boundz)
    boundam = np.array(boundam)
    boundhp = np.array(boundhp)
    boundrho = np.array(boundrho)
    boundvx = np.array(boundvx)
    boundvy = np.array(boundvy)
    boundvz = np.array(boundvz)
    boundvxdot = np.array(boundvxdot)
    boundvydot = np.array(boundvydot)
    boundvzdot = np.array(boundvzdot)
    boundu = np.array(boundu)
    boundudot = np.array(boundudot)
    boundgrpot = np.array(boundgrpot)
    boundmeanmolecular = np.array(boundmeanmolecular)
    boundcc = np.array(boundcc)
    bounddivv = np.array(bounddivv)
    bound_jrot=np.array(bound_jrot)
    bound_temp=np.array(bound_temp)
    bound_pressure=np.array(bound_pressure)

    # h1 = np.array(h1)
    # he3 = np.array(he3)
    # he4 = np.array(he4)
    # c12 = np.array(c12)
    # n14 = np.array(n14)
    # o16 = np.array(o16)
    # ne20 = np.array(ne20)
    # mg24 = np.array(mg24)
    
    try:
        bound_data = {
            'ntot':len(boundx),
            'x':boundx,
            'y':boundy,
            'z':boundz,
            'am':boundam,
            'hp':boundhp,
            'rho':boundrho,
            'vx':boundvx,
            'vy':boundvy,
            'vz':boundvz,
            'vxdot':boundvxdot,
            'vydot':boundvydot,
            'vzdot':boundvzdot,
            'u':boundu,
            'udot':boundudot,
            'grpot':boundgrpot,
            'meanmolecular':boundmeanmolecular,
            'cc':boundcc,
            'divv':bounddivv,
            # 'aa':boundaa,
            # 'bb':boundbb,
            # 'dd':bounddd
            'cc1val':cc1val,
            'cc2val':cc2val,
            'n1':n1,
            'n2':n2
        }
    except:
        bound_data = {
            'ntot':len(boundx),
            'x':boundx,
            'y':boundy,
            'z':boundz,
            'am':boundam,
            'hp':boundhp,
            'rho':boundrho,
            'vx':boundvx,
            'vy':boundvy,
            'vz':boundvz,
            'vxdot':boundvxdot,
            'vydot':boundvydot,
            'vzdot':boundvzdot,
            'u':boundu,
            'udot':boundudot,
            'grpot':boundgrpot,
            'meanmolecular':boundmeanmolecular,
            'cc':boundcc,
            'divv':bounddivv
            # 'aa':boundaa,
            # 'bb':boundbb,
            # 'dd':bounddd
           # 'cc1val':cc1val,
           # 'cc2val':cc2val,
           # 'n1':n1,
           # 'n2':n2
        }        

    print("Number of bound particles to the star =", len(boundrho))
    # print(boundrho)

    bound_composition = [boundh1, boundhe3, boundhe4, boundc12, boundn14, boundo16, boundne20, boundmg24]

    return bound_data,bound_composition,bound_jrot,bound_temp,bound_pressure
