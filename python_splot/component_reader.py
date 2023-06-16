"""
This routine quickly reads in the a file that contains
the x, y, z, and component data of each particle from
a collision simulation. It then saves every particle that
belongs to the star (with component 1) to a new file to
compare to the out*.sph file for continued analysis.

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
        FNAME = f'comp000{component_file_number}.sph'
    elif component_file_number < 100:
        FNAME = f'comp00{component_file_number}.sph'
    elif component_file_number < 1000:
        FNAME = f'comp0{component_file_number}.sph'
    else:
        FNAME = f'comp{component_file_number}.sph'

    with open(FNAME, "r") as f:

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
    
    # with open("starbound.sph", "w") as file:
    #     for i in range(len(component)):
    #         file.write(f"{x[i]}     {y[i]}     {z[i]}\n")
    
    # file.close()

    return component
