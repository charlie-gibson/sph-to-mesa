def write_bound_particles(nnit, header_data, bound_data, starnum):

    """
    nnit = input number
    header_data = header data from out****.sph
    bound_data = SPH data bound to a specific component
    starnum = the component number that we want the SPH particles to be bound to
    
    This code uses just the particle data bound to a given
    component (star 1, star 2, etc.) and writes a binary file
    that can be used to intialize a new sph simulation

    Support currently only exists for finding the particles of a single
    bound state. I will likely make an update for getting the bound
    particles of multiple bound states (e.g. starnum = 1 and 2)

    Charles Gibson
    Department of Physics
    Allegheny College
    01/01/2025
    """

    import numpy as np
    import struct

    # gets particle data for only bound particles
    data = {
        'x':np.array(bound_data['x']),
        'y':np.array(bound_data['y']),
        'z':np.array(bound_data['z']),
        'am':np.array(bound_data['am']),
        'hp':np.array(bound_data['hp']),
        'rho':np.array(bound_data['rho']),
        'vx':np.array(bound_data['vx']),
        'vy':np.array(bound_data['vy']),
        'vz':np.array(bound_data['vz']),
        'vxdot':np.array(bound_data['vxdot']),
        'vydot':np.array(bound_data['vydot']),
        'vzdot':np.array(bound_data['vzdot']),
        'u':np.array(bound_data['u']),
        'udot':np.array(bound_data['udot']),
        'grpot':np.array(bound_data['grpot']),
        'meanmolecular':np.array(bound_data['meanmolecular']),
        'cc':np.array(bound_data['cc']),
        'divv':np.array(bound_data['divv'])# ,
#        'aa':np.array(data['aa'])[inds],
#        'bb':np.array(data['bb'])[inds],
#        'dd':np.array(data['dd'])[inds],
#        'junk':np.full(new_ntot, 2.97079411e-312)
}

    new_ntot = len(data['x'])
    print('NEW NTOT:',new_ntot)

    # overwrites ntot with the new ntot
    header_data[0]=new_ntot
#    headers.append(0)

    header_data = tuple(header_data)

    # specifies how each value is stored in the header data
    header_dtype = np.dtype([
        ('ntot', np.int32), ('nnopt', np.int32), ('hco', np.float64), ('hfloor', np.float64),
        ('sep0', np.float64), ('tf', np.float64), ('dtout', np.float64), ('nout', np.int32),
        ('nit', np.int32), ('t', np.float64), ('nav', np.int32), ('alpha', np.float64),
        ('beta', np.float64), ('tjumpahead', np.float64), ('ngr', np.int32), ('nrelax', np.int32),
        ('trelax', np.float64), ('dt', np.float64), ('omega2', np.float64), ('ncooling', np.int32),
        ('erad', np.float64), ('ndisplace', np.int32), ('displacex', np.float64), ('displacey', np.float64),
        ('displacez', np.float64) #, ('junk', np.float64)
    ])

    # specifies how each value is stored in the particle data
    particle_dtype = np.dtype([
        ('x', np.float64), ('y', np.float64), ('z', np.float64), ('am', np.float64),
        ('hp', np.float64), ('rho', np.float64), ('vx', np.float64), ('vy', np.float64),
        ('vz', np.float64), ('vxdot', np.float64), ('vydot', np.float64), ('vzdot', np.float64),
        ('u', np.float64), ('udot', np.float64), ('grpot', np.float64), ('meanmolecular', np.float64),
        ('cc', np.int32), ('divv', np.float64) #,
        #('aa', np.float64), ('bb', np.float64), ('dd', np.float64),
#        ('junk', np.float64)
    ], align=False)

    # debug statement
#    print([type(value) for value in header_data])

    body_data = list(zip(
        data['x'], data['y'], data['z'], data['am'], data['hp'], data['rho'],
        data['vx'], data['vy'], data['vz'], data['vxdot'], data['vydot'], data['vzdot'],
        data['u'], data['udot'], data['grpot'], data['meanmolecular'], data['cc'], data['divv'] #,
        #data['aa'], data['bb'], data['dd'],
#        data['junk']
    ))

    # Prepare header and particle data as byte arrays
    header_array = np.array(tuple(header_data), dtype=header_dtype)
    particle_array = np.zeros(new_ntot, dtype=particle_dtype)
    for key, values in data.items():
        particle_array[key] = values

#    print("Particle data:", particle_array[:5])
#    raise SystemExit

    # Write the data to a binary file
    with open(f"bound{nnit:04d}_{starnum}.sph", "wb") as f:

        # Fortran innately writes marker data at the beginning and end of each line
        # Python does not, so we need to manually add this data instead

        # gets the number of bytes of the header data,
        # particle data, and the new ntot
        header_bytes = header_array.tobytes()
        particle_bytes = particle_array.tobytes()
        ntot_bytes = struct.pack('i', new_ntot)

        # Writes header data
        #header_bytes = struct.pack(header_dtype, *header_data)
        f.write(struct.pack('i', len(header_bytes)))  # Marker (size)
        f.write(header_bytes)
        f.write(struct.pack('i', len(header_bytes)))  # Marker (size)
    
        # Write particle records
        for i in range(len(particle_array)):
            particle_bytes_i = particle_array[i].tobytes()# struct.pack(particle_dtype[i], *particle_array[i])
            f.write(struct.pack('i', len(particle_bytes_i)))  # Marker (size)
            f.write(particle_bytes_i)
            f.write(struct.pack('i', len(particle_bytes_i)))  # Marker (size)

        # Write footer
        footer_bytes = struct.pack('i', new_ntot)
        f.write(struct.pack('i', len(footer_bytes)))  # Marker (size)
        f.write(footer_bytes)
        f.write(struct.pack('i', len(footer_bytes)))  # Marker (size)
