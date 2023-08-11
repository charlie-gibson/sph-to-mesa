"""
This code reads in the data of the particles from each zone (shell) of the star
from the profile exported by MESA. You will need to edit the name of the profile for
your own trials.

This file should be run on a relaxed model of the star, not a perturbed star

Charles Gibson
Allegheny College
Department of Physics
05/26/2023
"""

def composition_reader(profile_num):

    import numpy as np

    # profile15 is the default value
    # change if you are using a different profile
    with open(f"profile{profile_num}.data", "r") as f:
        # mass fraction
        q = []

        # composition
        h1 = []
        he3 = []
        he4 = []
        c12 = []
        n14 = []
        o16 = []
        ne20 = []
        mg24 = []

        # indexing value
        n = 0

        for line in f:
            # skips the first 6 lines as they don't contain pertinent information for us
            if n < 6:
                n += 1
            if n == 6:
                # determines the header data and finds the correct column to pull data from to append
                headers = line.split()
                q_index = headers.index("q")
                h1_index = headers.index("h1")
                he3_index = headers.index("he3")
                he4_index = headers.index("he4")
                c12_index = headers.index("c12")
                n14_index = headers.index("n14")
                o16_index = headers.index("o16")
                ne20_index = headers.index("ne20")
                mg24_index = headers.index("mg24")
                n+=1
            elif n > 6:
                values = line.split()
                # finds the mass fraction value
                qval = float(values[q_index])
                q.append(qval)

                # appends the value for each element to the corresponding position in the array
                h1val = float(values[h1_index])
                he3val = float(values[he3_index])
                he4val = float(values[he4_index])
                c12val = float(values[c12_index])
                n14val = float(values[n14_index])
                o16val = float(values[o16_index])
                ne20val = float(values[ne20_index])
                mg24val = float(values[mg24_index])

                h1.append(h1val)
                he3.append(he3val)
                he4.append(he4val)
                c12.append(c12val)
                n14.append(n14val)
                o16.append(o16val)
                ne20.append(ne20val)
                mg24.append(mg24val)

    f.close()

    q = np.array(q) # converts q (a list) to a numpy array for easier calculations in further codes

    elements = [h1, he3, he4, c12, n14, o16, ne20, mg24]

    return q, elements
