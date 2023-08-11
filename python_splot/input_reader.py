"""
This code reads in the sph.input file used for a StarSmasher
run and saves all of the data into a dictionary

Charles Gibson
Allegheny College
Department of Physucs
07/06/2023
"""

def sph_input_reader():

    # creates empty dictionary
    sph_data = {}
    
    # reads in file
    with open(sph.input, 'r') as f:

        n = 0
        
        for line in f:
            if n == 0 or n == 1:
                n+=1
            else:
                vals = line.split()
                var = vals[0]
                val = vals[2]
                try:
                    sph_data[var] = int(val)
                except try:
                    sph_data[var] = float(val)
                except:
                    sph_data[var] = val

    return sph_data
                
