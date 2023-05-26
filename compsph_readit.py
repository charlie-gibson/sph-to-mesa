"""
This code reads in all the data from the composition.sph file and saves it to
a dictionary to be read in by bestfit.py

By default, there are 8 values of composition on each line. You will need to
edit the function to alter the proper composition values

Charles Gibson
Allegheny College
Department of Physics
"""

def compsph_readit():

    # intializes element lists
    h1 = []
    he3 = []
    he4 = []
    c12 = []
    n14 = []
    o16 = []
    ne20 = []
    mg24 = []

    # reads composition.sph
    with open("composition.sph", "r") as f:

        # reads in each line, appending the correct composition data to the corresponding list
        for line in f:
            values = line.split()
            h1.append(float(values[0]))
            he3.append(float(values[1]))
            he4.append(float(values[2]))
            c12.append(float(values[3]))
            n14.append(float(values[4]))
            o16.append(float(values[5]))
            ne20.append(float(values[6]))
            mg24.append(float(values[7]))
    f.close()

    # creates a list for all data to be exported
    composition_data = [h1, he3, he4, c12, n14, o16, ne20, mg24]

    return composition_data
