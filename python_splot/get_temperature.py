"""
This code uses the density internal enregy outputted from StarSmasher
to determine the temperature of the particles.

It was originally written by Dr. Jamie Lombardi and Scott Fleming
in Fortran to be used in the src code of a StarSmasher run.
This code is just a translation from Fortran to Python.

The original commented heading is as follows:
---------------------------------------------------------------------
    subroutine to solve 4th order equations to determine the temperature x3
    for an equation of state with both ideal gas and radiation pressure.
    written by scott fleming 10/04/02 and james lombardi 2002-2003

    the fourth order equation comes from u_gas+ u_rad = u, with
    u_gas proportional to t and u_rad proportional to t^4

    in general, we can transform a 4th order equation to x^4+px^2+qx+r=0
    (see pages 57-58 of stillwell's "mathematics and its history" text)
    but we fortunately don't even have an x^2 term (that is, p=0).

    follow stillwell, we can transform this into a cubic equation:
    first solve for y by using the fact that b^2-4ac=0
    equation is then:  y^3=ry+q^2/8
    using the solution of cubic equations found in stillwell page 55:
"""

def get_temperature(q, r):

    import numpy as np

    k = 0.125 * q**2
    kh = 0.5*k

    if kh**2 - (r/3)**3 < 0:
        print(k, r, kh**2 - (r/3)**3)
        print("bad input: maybe imaginary results?")
        raise SystemExit

    piece1 = kh + (kh**2 - (r/3)**3)**0.5
    piece2 = (r/3)**3/piece1

    y1 = piece1**(1/3)

    y2 = -abs(piece2)**(1/3)
    yy = y1 + y2

    aa = 2*yy
    b = -q

    b2 = aa**0.5
    c2 = 0.5 * b/b2 + yy

    x3 = -2 * c2 / (b2 + (b2**2 - 4)**(0.5))

    if piece1 < 0:
        print(f"piece1 < 0", k, r, piece1, piece2)

    if b2 == 0:
        x3 = (-r - q*(-r-q*(-r)**0.25)**0.25)**0.25

    if piece2 >= 0:
        x3 = -(r + (r/q)**4) / q

    return x3
