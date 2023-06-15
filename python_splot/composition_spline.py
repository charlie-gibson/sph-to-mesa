"""
This code interpolates points between data to get a function of
mass fraction to the composition. It will yield the same number
of plots as there are elements kept track of in the star.

By default, there will be 8 elements in the order of:
H1
He3
He4
C12
N14
O16
Ne20
Mg24

If you do not have all of these outputted, then you will need to edit the code
to account for what elements are present.

Charles Gibson
Allegheny College
Department of Physics
05/26/2023
"""

def comp_spline(q, elements):

    import numpy as np
    from scipy.interpolate import CubicSpline
    import matplotlib.pyplot as plt

    h1 = elements[0]
    he3 = elements[1]
    he4 = elements[2]
    c12 = elements[3]
    n14 = elements[4]
    o16 = elements[5]
    ne20 = elements[6]
    mg24 = elements[7]

    q = q[::-1]
    h1.reverse()
    he3.reverse()
    he4.reverse()
    c12.reverse()
    n14.reverse()
    o16.reverse()
    ne20.reverse()
    mg24.reverse()

    # creates a cubic spline of each element vs the mass fraction to interpolate data
    h1fit = CubicSpline(q, h1)
    he3fit = CubicSpline(q, he3)
    he4fit = CubicSpline(q, he4)
    c12fit = CubicSpline(q, c12)
    n14fit = CubicSpline(q, n14)
    o16fit = CubicSpline(q, o16)
    ne20fit = CubicSpline(q, ne20)
    mg24fit = CubicSpline(q, mg24)

    q = q[::-1]
    h1.reverse()
    he3.reverse()
    he4.reverse()
    c12.reverse()
    n14.reverse()
    o16.reverse()
    ne20.reverse()
    mg24.reverse()

    comp_splines = {
        'h1':h1fit,
        'he3':he3fit,
        'he4':he4fit,
        'c12':c12fit,
        'n14':n14fit,
        'o16':o16fit,
        'ne20':ne20fit,
        'mg24':mg24fit
    }

    x = np.linspace(np.min(q), np.max(q), len(q))

    plt.plot(x, h1fit(x))
    plt.title("Hydrogen Bestfit", fontsize=9)
    plt.xlabel("$q$")
    plt.ylabel("$H1$")
    plt.show()

    return comp_splines
