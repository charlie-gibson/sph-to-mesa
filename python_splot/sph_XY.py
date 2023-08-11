"""
This code quickly finds the star's total X, Y, and Z fractions
to provide as a check after the composition file has been generated

Charles Gibson
08/11/2023
Department of Physics
Allegheny College
Center for Interdisciplinary Expolration and Research in Astrophysics (CIERA)
Northwestern University
"""

def xy_reader():
    import numpy as np
    import matplotlib.pyplot as plt

    with open('composition.dat','r') as f:
        n=0
        xq=[]
        h1=[]
        he4=[]
        for line in f:
            if n==0:
                n+=1
            else:
                cols=line.split()
                xq.append(float(cols[0]))
                h1.append(float(cols[1]))
                he4.append(float(cols[3]))

                weightedX=0
                weightedY=0

    for i in range(len(xq)-1):
        xq_frac=xq[i+1]-xq[i]
        weightedX+=h1[i]*xq_frac
        weightedY+=he4[i]*xq_frac

    Z=1-weightedX-weightedY

    print(f'X: {weightedX}\nY: {weightedY}\nZ: {Z}')
