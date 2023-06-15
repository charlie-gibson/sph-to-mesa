import numpy as np
import matplotlib.pyplot as plt

def read_eos():
    eosfile = "sph.eos_X0.30to0.75step0.05"
    print(f"about to read EOS file {eosfile}")

    with open(eosfile, "r") as f:
        numx, xtablefirst, xtablelast, stepx, *_ = f.readline().split()
        numx = int(numx)
        xtablefirst = float(xtablefirst)
        xtablelast = float(xtablelast)
        stepx = float(stepx)
        print(f"There are {numx} hydrogen abundances, ranging from {xtablefirst} to {xtablelast}, in steps of {stepx}")

        for ix in range(numx):
            xxx, *_ = f.readline().split()
            yyy, *_ = f.readline().split()
            zzz, *_ = f.readline().split()
            abar, *_ = f.readline().split()
            zbar, *_ = f.readline().split()
            xxx = float(xxx)
            yyy = float(yyy)
            zzz = float(zzz)
            abar = float(abar)
            zbar = float(zbar)

            eosmu = 1.67262158e-24 / (2 * xxx + 0.75 * yyy + 0.5 * zzz)

            print('x,y,z=', xxx, yyy, zzz)
            print('abar,zbar=',abar,zbar)
            print('EOS mu=', eosmu, 1.67262158e-24 *abar/(1+zbar))

            if abs(xxx + yyy + zzz - 1.0) > 1e-16:
                print('one minus 1=', xxx + yyy + zzz - 1.0)
                raise SystemExit

            numrho, rhotablefirst, rhotablelast, steprho = map(float, f.readline().split()[:4])
            numu, utablefirst, utablelast, stepu = map(float, f.readline().split()[:4])
            numrho = int(numrho)
            numu = int(numu)

            print('numrho, rhotablefirst, rhotablelast, steprho=', numrho, rhotablefirst, rhotablelast, steprho)
            print('numu, utablefirst, utablelast, stepu=',numu, utablefirst, utablelast, stepu)

            for i in range(2):
                f.readline().split()

            if (ix == 0):
                rhotable = np.empty(numrho, dtype = float)
                utable = np.empty(numu, dtype = float)
                eostable = np.empty((numu,numrho,numx,3), dtype = float)

            for irho in range(numrho):
                for iu in range(numu):
                    inputs = list(map(float, f.readline().split()))
                    rhotable[irho] = inputs[0]
                    utable[iu] = inputs[1]
                    eostable[iu][irho][ix][:] = inputs[2:5]
                    eostable[iu][irho][ix][0] = 10.0 ** eostable[iu][irho][ix][0] # record temperature instead of log temperature
                    eostable[iu][irho][ix][1] *= 1.67262158e-24 # record mean molecular mass in grams
                #ttable = eostable[irho,:,ix,0]
                #print('ttable=',ttable)
                #print('utable=',utable)
                #raise SystemExit
                #if irho % 40 == 0 and ix == numx-1:
                #    plt.plot(eostable[irho,:,ix,0], utable, label=f'log rho={rhotable[irho]}')
            #if ix == numx-1:
            #    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            #    plt.xlabel("log T")
            #    plt.ylabel("log u")
            #    plt.title(f"X={xxx}, Y={yyy}, Z={zzz}")
            #    plt.tight_layout()
            #    plt.show()  

            print('done reading eos table')
            print('numrho=', numrho, 'numu=', numu)

            if rhotablefirst != rhotable[0]:
                print('rhotablefirst mismatch', rhotablefirst, rhotable[0], rhotablefirst - rhotable[0])
                raise SystemExit

            if utablefirst != utable[0]:
                print('utablefirst mismatch')
                raise SystemExit

            if rhotablelast != rhotable[numrho-1]:
                print('rhotablelast mismatch')
                raise SystemExit

            if utablelast != utable[numu-1]:
                print('utablelast mismatch')
                raise SystemExit

            steprhotest = (rhotable[numrho-1] - rhotable[0]) / (numrho - 1)
            steputest = (utable[numu-1] - utable[0]) / (numu - 1)

            if abs(steprhotest - steprho) > 1e-9:
                print('steprho problem')
                raise SystemExit

            if abs(steputest - stepu) > 1e-9:
                print('stepu problem')
                raise SystemExit
    
    eos_table = {
        'rhotablefirst':rhotablefirst,
        'steprho':steprho,
        'steprho':steprho,
        'utablefirst':utablefirst,
        'stepu':stepu,
        'xtablefirst':xtablefirst,
        'stepx':stepx,
        'numx':numx,
        'numrho':numrho,
        'numu':numu,
        'eostable':eostable
    }
    
    return eos_table

def useeostable(eos_data, ucgs, rhocgs, xxx, which):
    # ucgs = internal energy per unit mass in cgs units
    # rhocgs = density in cgs units
    # xxx = hydrogen mass fraction X
    #     which=0 returns temperature
    #     which=1 returns mean molecular mass mu
    #     which=2 returns pressure

    #     utable[iu]=utable[0]+iu*stepu
    #     so, iu= (utable[iu]-utable[0])/stepu

    # uses the data from the read_eos() function
    # it imports this data, so that if the routine is used many times, the code
    # will not need to read the data the whole time
    rhotablefirst = eos_data['rhotablefirst']
    steprho = eos_data['steprho']
    utablefirst = eos_data['utablefirst']
    stepu = eos_data['stepu']
    xtablefirst = eos_data['xtablefirst']
    stepx = eos_data['stepx']
    numx = eos_data['numx']
    numrho = eos_data['numrho']
    numu = eos_data['numu']
    eostable = eos_data['eostable']
    # print(eostable.shape, "eostable test print")

    # finds the log of the density and internal energy to be used in the table
    log10rho = np.log10(rhocgs)
    log10u = np.log10(ucgs)

    # rounds density and internal energies to the nearest step in the table
    irho = int((log10rho - rhotablefirst) / steprho)
    iu = int((log10u - utablefirst) / stepu)
    
    #xxx = (1.67262158e-24 / particlemu + 0.25 * zzz - 0.75) / 1.25
    ix = min(int((xxx - xtablefirst) / stepx), numx - 2)
    
    xlow = xxx - (xtablefirst + ix * stepx)
    xhigh = xtablefirst + (ix + 1) * stepx - xxx

    # if the x value is too small or large for the table
    if (ix < 0 or ix > numx - 1) and numx > 1:
        print("Expand X range covered by EOS file")
        if ix < 0:
            ix = 0
        elif ix == numx - 1:
            ix = numx - 2
        else:
            raise SystemExit
    
    if 0 <= irho <= numrho - 2 and 0 <= iu <= numu - 2:
        rholow = log10rho - (rhotablefirst + irho * steprho)
        rhohigh = rhotablefirst + (irho+1) * steprho - log10rho
        
        # print("iu, numu =", iu, numu)
        if iu>=0 and iu<=numu-2:

            ulow = log10u - (utablefirst + iu * stepu)
            uhigh = utablefirst + (iu+1) * stepu - log10u
            
            f00 = rholow * ulow
            f10 = rhohigh * ulow
            f01 = rholow * uhigh
            f11 = rhohigh * uhigh

            if numx > 1:
                # print("iu, numu, irho, ix, which", iu, numu, irho, ix, which)
                useeostable0 = f00 * eostable[iu + 1, irho + 1, ix, which] \
                    +          f10 * eostable[iu + 1, irho, ix, which] \
                    +          f01 * eostable[iu, irho + 1, ix, which] \
                    +          f11 * eostable[iu, irho, ix, which]
                # useeostable0 = eostable[iu, irho, ix, which]
                useeostable1 = f00 * eostable[iu + 1][irho + 1][ix + 1][which] \
                    +          f10 * eostable[iu + 1][irho][ix + 1][which] \
                    +          f01 * eostable[iu][irho + 1][ix + 1][which] \
                    +          f11 * eostable[iu][irho][ix + 1][which]
                useeostable = (xlow * useeostable1 + xhigh * useeostable0) / (steprho * stepu * stepx)
            else:
                useeostable = (f00 * eostable[iu + 1][irho + 1][0][which] \
                    +          f10 * eostable[iu + 1][irho][0][which] \
                    +          f01 * eostable[iu][irho + 1][0][which] \
                    +          f11 * eostable[iu][irho][0][which] \
                ) / (steprho * stepu)
        elif iu <= 1:
            #     First we get the value at the edge of the table with smallest u 
            if numx > 1:
                useeostable0 = rholow * eostable[0][irho + 1][ix][which] \
                    +         rhohigh * eostable[0][irho][ix][which] 
                useeostable1 = rholow * eostable[0][irho + 1][ix + 1][which] \
                    +         rhohigh * eostable[0][irho][ix + 1][which]
                useeostable = (xlow * useeostable1 + xhigh * useeostable0) / (steprho * stepx)
            else:
                useeostable = (rholow * eostable[0][irho + 1][0][which] + rhohigh * eostable[0][irho][0][which]) / steprho
            
            if which != 1:
                #     Reminder: which=1 returns mean molecular weight from EOS table
                #     so if which != 1 the code is supposed to calclate pressure or temperature.
                #     We are at very low specific internal energy u, where the pressure
                #     or temperature should be nearly proportional to rho*u (at fixed
                #     composition).  We use this to extrapolate to smaller u:
                useeostable = ucgs / 10 ** utablefirst * useeostable
        elif iu >= numu-1:
            if which == 2:
                #     We are at very high specific internal energy u, where the pressure
                #     nearly proportional to u (at fixed rho and composition).     
                #     Use linear interpolation between the two cartesian
                #     grid points (irho,numu), (irho+1,numu) and in addition extrapolate
                #     to larger u.
                if numx > 1:
                    useeostable0 = rholow * eostable[numu-1][irho + 1][ix][which] \
                        +          rhohigh * eostable[numu-1][irho][ix][which]
                    useeostable1 = rholow * eostable[numu-1][irho + 1][ix + 1][which] \
                        +          rhohigh * eostable[numu-1][irho][ix + 1][which]
                    useeostable = ucgs / 10 ** (utablefirst + (numu - 1) * stepu) * ( \
                        xlow * useeostable1 + xhigh * useeostable0 \
                    ) / (steprho * stepx)
                else:
                    useeostable = (rholow * eostable[numu-1][irho + 1][0][which] \
                        +         rhohigh * eostable[numu-1][irho][0][which] \
                    ) / steprho
            elif which == 0:
                #     We are at very high specific internal energy u, where the temperature
                #     is nearly proportional to u^(1/4) (at fixed rho and composition)     .         
                #     Use linear interpolation between the two cartesian
                #     grid points (irho,numu), (irho+1,numu) and in addition extrapolate
                #     to larger u.

                if numx > 1:
                    useeostable0 = rholow * eostable[numu-1][irho + 1][ix][which] \
                        +         rhohigh * eostable[numu-1][irho][ix][which]
                    useeostable1 = rholow * eostable[numu-1][irho + 1][ix + 1][which] \
                        +         rhohigh * eostable[numu-1][irho][ix + 1][which]
                    useeostable = ((ucgs / 10 ** (utablefirst + (numu - 1) * stepu)) ** 0.25 \
                        * (xlow * useeostable1 + xhigh * useeostable0) \
                    ) / (steprho * stepx)
                else:
                    useeostable = (rholow * eostable[numu-1][irho + 1][0][which] \
                        +         rhohigh * eostable[numu-1][irho][0][which] \
                    ) / steprho
            else:
                if numx > 1:
                    useeostable0 = rholow * eostable[numu-1][irho + 1][ix][which] \
                        +         rhohigh * eostable[numu-1][irho][ix][which]
                    useeostable1 = rholow * eostable[numu-1][irho + 1][ix + 1][which] \
                        +         rhohigh * eostable[numu-1][irho][ix + 1][which]
                    useeostable = (xlow * useeostable1 + xhigh * useeostable0) / (steprho * stepx)
                else:
                    useeostable = (rholow * eostable[numu-1][irho + 1][0][which] \
                        +         rhohigh * eostable[numu-1][irho][0][which] \
                    ) / steprho
    else:
        #     at extreme densities we will use ideal gas + radiation pressure
        if irho < 0:
            irho = 0
        if irho > numrho-1:
            irho = numrho-1
            
        if iu < 0:
            if numx > 1:
                meanmu0 = eostable[0][irho][ix][1]
                meanmu1 = eostable[0][irho][ix + 1][1]
                meanmu = (xlow * meanmu1 + xhigh * meanmu0) / stepx
            else:
                meanmu = eostable[0][irho][0][1]
        elif iu >= numu-1:
            if numx > 1:
                meanmu0 = eostable[numu-1][irho][ix][1]
                meanmu1 = eostable[numu-1][irho][ix + 1][1]
                meanmu = (xlow * meanmu1 + xhigh * meanmu0) / stepx
            else:
                meanmu = eostable[numu-1][irho][0][1]
        else:
            ulow = log10u - (utablefirst + iu * stepu)
            uhigh = utablefirst + (iu + 1) * stepu - log10u
            if numx > 1:
                meanmu0 = ulow * eostable[iu + 1][irho][ix][1] + uhigh * eostable[iu][irho][ix][1]
                meanmu1 = ulow * eostable[iu + 1][irho][ix + 1][1] + uhigh * eostable[iu][irho][ix + 1][1]
                meanmu = (xlow * meanmu1 + xhigh * meanmu0) / (stepu * stepx)
            else:
                meanmu = (ulow * eostable[iu + 1][irho][0][1] + uhigh * eostable[iu][irho][0][1]) / stepu
        
        if which == 1:
            useeostable = meanmu
        else:
            # We will get the temperature from the ideal gas + radiation pressure analytic EOS. Note:
            # ionization energy is neglected, which could be bad at low temperatures.

            # this analytical method can be found in get_temperature.py
            print("WARNING: OUTSIDE THE TABLE BOUNDARY ON COLD END. DO NOT TRUST THE T or P")

            # helpful constants
            boltz=1.380658e-16 # erg/kelvin
            crad=2.997924580e10 # cm/sec  NOTE: card has a different meaning in MESA
            planck=6.6260755e-27 # gram cm^2/sec
            crad2=crad**2 # cm^2/s^2
            pi=3.1415926535897932384626
            sigma=pi**2*boltz*(boltz*2*pi/planck)**3/60/crad2 #cgs
            arad=4.0*sigma/crad #cgs
            qconst=1.5*boltz/arad #cgs
            q=qconst*rhocgs/meanmu
            r=-ucgs*rhocgs/arad
            k = 0.125 * q ** 2
            kh = 0.5 * k
            piece1 = kh + np.sqrt(kh ** 2 - (r / 3) ** 3)
            piece2 = (r / 3) ** 3 / piece1
            y1 = piece1 ** (1 / 3)
            y2 = -abs(piece2) ** (1 / 3)
            yy = y1 + y2
            aa = 2 * yy
            b = -q
            b2 = np.sqrt(aa)
            if b2 == 0:
                temperature=(-r -q*( -r-q*(-r)**0.25 )**0.25)**0.25
                print('hi', meanmu)
            else:
                c2 = 0.5 * b / b2 + yy
                temperature = -2 * c2 / (b2 + np.sqrt(b2 ** 2 - 4 * c2))
            if which == 0:
                useeostable = temperature
            elif which ==2:
                pgas=rhocgs*boltz*temperature/meanmu
                prad=arad*temperature**4/3
                beta1=pgas/(pgas+prad)
                gam1=(32-24*beta1-3*beta1**2) / (24-21*beta1)
                
                if gam1 <= 0.999*4/3 or gam1 >= 1.001*5/3:
                    print('error gam1=',gam1)
                    print(beta1,pgas,prad,temperature,-ucgs*rhocgs/arad)
                    raise SystemExit
                useeostable=pgas+prad
    
    return useeostable

# Let's test the function
# logrho = -14
# xxx = 0.3
# print(f"X={xxx}")
# print("logrho logu T mu P")
# for logu in np.arange(9.8, 9.9, 0.01):
#     print(logrho, logu, useeostable(ucgs=10**logu, rhocgs=10**logrho, xxx=xxx, which=0), \
#                         useeostable(ucgs=10**logu, rhocgs=10**logrho, xxx=xxx, which=1), \
#                         useeostable(ucgs=10**logu, rhocgs=10**logrho, xxx=xxx, which=2))
