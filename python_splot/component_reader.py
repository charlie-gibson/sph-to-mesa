"""
This routine reads in the the comp****.sph file that contains
the x, y, z, and bound component data of each particle from
a collision simulation. It then returns the bound state of each
particle, where
1 = bound to star 1
2 = bound to star 2
3 = bound to star 3
4 = unbound

Charles Gibson
Allegheny College
Department of Physics
"""

def component_reader(component_file_number:str):

    x = []
    y = []
    z = []
    component = []

    # used to ensure correct naming of the file
    if component_file_number < 10:
        fname = f'comp000{component_file_number}.sph'
    elif component_file_number < 100:
        fname = f'comp00{component_file_number}.sph'
    elif component_file_number < 1000:
        fname = f'comp0{component_file_number}.sph'
    else:
        fname = f'comp{component_file_number}.sph'

    with open(fname, "r") as f:

        n = 0

        for line in f:
            if n == 0:
                n += 1
            elif n > 0:
                values = line.split()
                xval = float(values[0])
                yval = float(values[1])
                zval = float(values[2])
                componentval = int(values[3])
                component.append(componentval)
                x.append(xval)
                y.append(yval)
                z.append(zval)
    f.close()

    return component
