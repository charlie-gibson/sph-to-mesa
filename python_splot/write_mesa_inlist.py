"""
Writes a MESA inlist from SPH data

Charles Gibson
Allegheny College
Department of Physics
12/04/2024
"""

def write_mesa_inlist(sph_data, comp_data, mode, starnum, hse = False):

    import numpy as np
    import os

    G     = 6.6743e-8
    munit = 1.9891e33
    runit = 6.9599e10

    # gets necessary star data
    am = sph_data['am']
    mass = np.sum(am)

    x = sph_data['x']
    y = sph_data['y']
    z = sph_data['z']

    # center of mass
    xc = np.sum(am * x) / mass
    yc = np.sum(am * y) / mass
    zc = np.sum(am * z) / mass
    
    # gets dynamical timescale of the star for an initial estimate
    # of entropy relax timescale

    r = np.sqrt((x - xc)**2 + (y - yc)**2 + (z - zc)**2)

    R = np.max(r) * runit
    print(R)
    M = mass * munit

    dyn_t = np.sqrt(R**3 / (G * M)) # in seconds
    dyn_t = dyn_t * 1/60 * 1//60 * 1/24 * 1/365

    dyn_t_str = f'{dyn_t / 1000:.2e}'.replace('e','d')

    print(f'DYNAMICAL TIMESCALE:  {dyn_t} yrs')

    h1 = np.array(comp_data[0])
    he4 = np.array(comp_data[2])

    m_h1 = np.sum(am * h1)
    m_he4 = np.sum(am * he4)

    z_tot = 1 - (m_h1 + m_he4) / mass

    ### writes inlist project with main parameters
    with open('inlist_project','w') as f:

        # writes star_job information
        f.write(f'&star_job\n')
        f.write(f'    create_pre_main_sequence_model = .false.\n\n')
        f.write(f'    save_model_when_terminate = .false.\n')
        f.write(f'    save_model_filename = \'{mass:.4f}Msun.mod\'\n\n')
        f.write(f'    pgstar_flag = .true.\n\n')

        ## input files from StarSmasher output
        f.write(f'  ! relaxing from StarSmasher output\n')
        f.write(f'  ! relax composition\n')
        f.write(f'    relax_initial_composition = .true.\n')
        if hse:
            f.write(f'    relax_composition_filename = \'composition{starnum}_hse.dat\'\n')
        else:
            f.write(f'    relax_composition_filename = \'composition{starnum}.dat\'\n')
        f.write(f'    num_steps_to_relax_composition = 100\n\n')
        f.write(f'  ! relax angular momentum\n')
        f.write(f'    relax_initial_angular_momentum = .true.\n')
        if hse:
            f.write(f'    relax_angular_momentum_filename = \'angular_momentum{starnum}_hse.dat\'\n')
        else:
            f.write(f'    relax_angular_momentum_filename = \'angular_momentum{starnum}.dat\'\n')
        f.write(f'    max_steps_to_relax_angular_momentum = 100\n\n')
        f.write(f'  ! relax entropy\n')
        f.write(f'    relax_initial_entropy = .true.\n')
        if hse:
            f.write(f'    relax_entropy_filename = \'entropy{starnum}_hse.dat\'\n')
        else:
            f.write(f'    relax_entropy_filename = \'entropy{starnum}.dat\'\n')
        f.write(f'  ! timescale used in Eq. (82) of Paxton et al. (2018)\n')
        ## writes the dynamical timescale based on the dynamical timescale:
        ##
        ## from Paxton et al. (2018)
        ## "The value τ should be chosen to be small
        ## enough that energy transport is negligible during the pseudoevolution. In practice, τ can be chosen to be orders of
        ## magnitude smaller than the dynamical timescale of the system."
        ##
        ## I arbitrarily choose τ to be the dynamical timescale dyn_t / 1000
        ## however this will not likely yield a reasonable model, requiring some changes to the value
        f.write(f'    timescale_for_relax_entropy = {dyn_t_str} ! change this value if there are convergence issues during relaxation\n')
        f.write(f'    max_steps_to_relax_entropy = 1000\n')
        if mode == 'S':
            pass
        else:
            f.write(f'    get_entropy_for_relax_from_eos = \'eos{mode}\'\n')
        f.write(f'/ ! end of star_job namelist\n\n')

        # writes eos information
        f.write(f'&eos\n')
        f.write(f'/ ! end of eos namelist\n\n')

        # writes kap information
        f.write(f'&kap\n')
        f.write(f'    use_Type2_opacities = .true.\n')
        f.write(f'    Zbase = {z_tot:.4f} ! may want to use Z in outer layers of star rather than total. Only really applies for he4 burning stars and older\n')
        f.write(f'/ ! end of kap namelist\n\n')

        # writes controls information
        f.write(f'&controls\n')
        f.write(f'    initial_mass = {mass:.4f} ! in Msun units\n')
        f.write(f'  ! when to stop\n')
        f.write(f'  ! when to stop\n')
        f.write(f'  ! stop when the star nears ZAMS (Lnuc/L > 0.99)\n')
        f.write(f'    Lnuc_div_L_zams_limit = 0.99d0\n')
        f.write(f'    stop_near_zams = .false.\n')
        f.write(f'  ! stop when the center mass fraction of h1 drops below this limit\n')
        f.write(f'    xa_central_lower_limit_species(1) = \'he4\' ! you can change this to other elements. for example h1, he3, c12, n14, o16, ne20, mg24\n')
        f.write(f'    xa_central_lower_limit(1) = 1d-5\n\n')
        f.write(f'  ! wind\n\n\n')
        f.write(f'  ! atmosphere\n\n\n')
        f.write(f'  ! rotation\n\n\n')
        f.write(f'  ! element diffusion\n\n\n')
        f.write(f'  ! mlt\n\n\n')
        f.write(f'  ! mixing\n\n\n')
        f.write(f'  ! timesteps\n')
        f.write(f'      min_timestep_limit = 1d-10\n\n')
        f.write(f'  ! mesh\n\n\n')
        f.write(f'  ! solver\n')
        f.write(f'    ! options for energy conservation (see MESA V, Section 3)\n')
        f.write(f'      energy_eqn_option = \'dedt\'\n')
        f.write(f'      use_gold_tolerances = .true.\n\n')
        f.write(f'  ! output\n')
        f.write(f'      history_interval = 1\n')
        f.write(f'      profile_interval = 10\n')
        f.write(f'      profile_model = 1 ! choose a profile that you are interested in using, perhaps for analysis or for another StarSmasher relaxation\n\n')
        f.write(f'/ ! end of controls namelist\n')

    ### writes inlist
    with open('inlist','w') as f:

        f.write('! This is the first inlist file that MESA reads when it starts.\n\n')
        f.write('! This file tells MESA to go look elsewhere for its configuration\n')
        f.write('! info. This makes changing between different inlists easier, by\n')
        f.write('! allowing you to easily change the name of the file that gets read.\n\n')

        # star job
        f.write('&star_job\n\n')
        f.write('    read_extra_star_job_inlist(1) = .true.\n')
        f.write('    extra_star_job_inlist_name(1) = \'inlist_project\'\n\n')
        f.write('/ ! end of star_job namelist\n\n\n')

        # eos
        f.write('&eos\n\n')
        f.write('    read_extra_eos_inlist(1) = .true.\n')
        f.write('    extra_eos_inlist_name(1) = \'inlist_project\'\n\n')
        f.write('/ ! end of eos namelist\n\n\n')

        # kap
        f.write('&kap\n\n')
        f.write('    read_extra_kap_inlist(1) = .true.\n')
        f.write('    extra_kap_inlist_name(1) = \'inlist_project\'\n\n')
        f.write('/ ! end of kap namelist\n\n\n')

        # controls
        f.write('&controls\n\n')
        f.write('    read_extra_controls_inlist(1) = .true.\n')
        f.write('    extra_controls_inlist_name(1) = \'inlist_project\'\n\n')
        f.write('/ ! end of controls namelist\n\n\n')

        #pgstar
        f.write('&pgstar\n\n')
        f.write('    read_extra_pgstar_inlist(1) = .true.\n')
        f.write('    extra_pgstar_inlist_name(1) = \'inlist_pgstar\'\n\n')
        f.write('/ ! end of pgstar namelist\n')

    with open('inlist_pgstar','w') as f:

        f.write('&pgstar\n')
        f.write('  ! see star/defaults/pgstar.defaults\n\n')
        f.write('  ! MESA uses PGPLOT for live plotting and gives the user a tremendous\n')
        f.write('  ! amount of control of the presentation of the information.\n\n')
        f.write('  ! show HR diagram\n')
        f.write('  ! this plots the history of L,Teff over many timesteps\n')
        f.write('    HR_win_flag = .true.\n\n')
        f.write('  ! set static plot bounds\n')
        f.write('    HR_logT_min = 3.5\n')
        f.write('    HR_logT_max = 4.6\n')
        f.write('    HR_logL_min = 2.0\n')
        f.write('    HR_logL_max = 6.0\n\n')
        f.write('  ! set window size (aspect_ratio = height/width)\n')
        f.write('    HR_win_width = 6\n')
        f.write('    HR_win_aspect_ratio = 1.0\n\n\n')
        f.write('  ! show temperature/density profile\n')
        f.write('  ! this plots the internal structure at single timestep\n')
        f.write('    TRho_Profile_win_flag = .true.\n\n')
        f.write('  ! add legend explaining colors\n')
        f.write('    show_TRho_Profile_legend = .true.\n\n')
        f.write('  ! display numerical info about the star\n')
        f.write('    show_TRho_Profile_text_info = .true.\n\n')
        f.write('  ! set window size (aspect_ratio = height/width)\n')
        f.write('    TRho_Profile_win_width = 8\n')
        f.write('    TRho_Profile_win_aspect_ratio = 0.75\n\n')
        f.write('/ ! end of pgstar namelist\n')

    # moves all MESA input files to a new directory called MESA_input_files
    
    os.system(f'mkdir MESA_input_files{starnum}')
    os.system(f'mv inlist inlist_project inlist_pgstar MESA_input_files{starnum}')
    if hse:
        os.system(f'cp -p angular_momentum{starnum}_hse.dat entropy{starnum}_hse.dat composition{starnum}_hse.dat MESA_input_files{starnum}')
    elif not hse:
        os.system(f'cp -p angular_momentum{starnum}.dat entropy{starnum}.dat composition{starnum}.dat MESA_input_files{starnum}')
