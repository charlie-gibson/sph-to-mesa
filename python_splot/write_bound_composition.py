def write_bound_composition(nnit, bound_composition, starnum):

    """
    Takes the bound particle data and writes it to a new
    bound sph.composition file

    parameters:
    nnit = correspond out****.sph file number
    bound_composition = the composition data bound to a given component
    starnum = the component number that the particles are bound to

    Charles Gibson
    Department of Physics
    Allegheny College
    01/01/2025
    """

    h1   = bound_composition[0]
    he3  = bound_composition[1]
    he4  = bound_composition[2]
    c12  = bound_composition[3]
    n14  = bound_composition[4]
    o16  = bound_composition[5]
    ne20 = bound_composition[6]
    mg24 = bound_composition[7]

    print('WRITING NEW SPH.COMPOSITION FILE')

    with open(f'sph.{nnit:04d}_bound_composition_{starnum}','w') as f:

        for h1val, he3val, he4val, c12val, n14val, o16val, ne20val, mg24val in zip(h1,
                                                                                he3,
                                                                                he4,
                                                                                c12,
                                                                                n14,
                                                                                o16,
                                                                                ne20,
                                                                                mg24):
            f.write(f'{h1val}     {he3val}     {he4val}     {c12val}     {n14val}     {o16val}     {ne20val}     {mg24val}\n')
