# sph-to-mesa description
This repository provides the code necessary to take an SPH model from StarSmasher and generate the necessary files to continue its evolution in MESA. It also provides a brief set of steps needed to use the data in MESA

Before you can use this code, you will need to make sure that the codes are set up for your machine.
  All codes found in python_splot should be in one location (you will not need multiple copies of these codes).
  Edit splot.py on lines 19 and 20 to include the python_splot directory location. splot.py should then be copied to whatever SPH directory you are working on.

  Ensure that the output in your StarSmasher (or other SPH code) is the same as here. For StarSmasher, in your simulation directory, go into src and look at output.f.
  If the values being read in by readit_collision.py are inconsistent, make sure that you make the necessary changes.

  Once you have made these adjustments, you should be ready to use the codes

Here is a description of the process to utilize the codes once they have been set up:
1. Ensure that the data from the original data mesa profile (i.e. profile15.data) contains all composition data for H1, He3, He4, C12, N14, O16, Ne20, and Mg24 along with the q (mass fraction) values
2. Relax your star using StarSmasher. See https://jalombar.github.io/starsmasher/ for more information
3. In your relaxation directory, run splot.py and choose option 2. This will create a file called composition.sph
5. Run the SPH simulation in a new directory and make sure that you know which output file(s) (out*.sph) you will be analyzing
6. Copy composition.sph to this new collision directory
    *** If you are using two stars, you will need to combine the two files into one with the composition.sph file for start1.u at the beginning and the composition.sph file for start2.u at the end
8. Run splot.py and use option 4. This will create a file called comp*.sph (where * is the file number and comp is short for component, not composition)
9. Run splot.py and use option 3. This will write 3 files to the directory: entropy.dat, composition.dat, and angular_momentum.dat
10. These three files can be used in the MESA as input files for the entropy, composition, and angular momentum

# Summary of Codes

## splot.py
* This code, originally created by Hans Vanderzyden and later modified by myself, acts as the main interface for the data analysis
* It's location should be in your working directory and it will import the other analysis files from their file locations
* There are currently two options for the code. Additional codes can be written for further analysis needs

The following scripts are for generating the necessary MESA files:

## readit_collision.py
* This code, originally written by Hans Vanderzyden and later modified by myself, reads in the binary data from the out*.sph file(s) and saves all the data to a dictionary to be exported
* This code is generalized to read in all particles except for BH particles (work may need to be done to ensure that this method of skipping the BH is robust)

## compsph_readit.py
* This code reads in composition data from the composition.sph file generated from the relaxed model
* The data is then saved as a list of lists which can be later read in in bestfit_total.py

## bestfit.py
* This code, translated from a Fortran routine by Dr. James Lombardi, uses the data from readit.py and bins all the data based on their densities
* Following the binning process, the script determines various quantities including mass (m), radius (r), pressure (P), temperature (T), specific internal energy (u), density (rho), specific angular momentum (j), and composition of each bin
* These values are all written to a new file called bestfit.sph

## composition_jpt_entropy_reader.py
* This code reads in bestfit.sph from bestfit_total.py and creates a spline of the data, specifically of u vs x, rho vs x, j vs x, and composition vs x 
* Then it reads in a dummy file, angular_momentum.dat and assigns reasonable mass fraction (q) values to each spline
* Finally, it writes the data points of each spline and the corresponding q value to three new files, entropy.dat, composition.dat, and angular_momentum.dat
* The file is formatted to have the first line have the number of zones in the new star (and the number of elements in composition.dat) followed by a new line for every zone, starting with the q value, followed by the appropriate value to be analyzed. angular_momentum.dat has the specific angular momentum (q   j), and composition.dat has the element abundances (q   H1   He3   He4   C12   N14   O16   Ne20   Mg24). The entropy can have several formats depending on the input. These include (q, rho, T), (q, P, T), or (q, rho, u).

The following scripts are for generating the correct composition files:

## read_orig_profile_comp.py
* This code reads in the composition of each element directly from the profile used in StarSmasher, saving each value to a list
* The code returns two lists. The first is the array of all q values that correspond with the shell, and the second is a list of lists that contains the composition for each corresponding q value

## composition_spline.py
* This code uses the data from read_orig_comp.py and creates a spline of the data for each element
* This allows data to be interpolated between each data point corresponding to the original zones

## composition_fit.py
* This code uses composition_spline.py and the data from readit.py to now find the composition of each particle of the SPH model using a similar binning approach to what bestfit_total.py does
* It then writes all the data to composition.sph. Every line contains new composition data for each particle in the same order that each particle is outputted in the SPH simulation

The following codes are used to determine the temperature for each particle (as it's not included in the output of the SPH)

## get_temperature.py
* This code uses an analytic, fourth degree solution to calculate the temperature based on the temperature and internal energy of a particle
* This is the same method used in the analytic EOS used in the StarSmasher

## eos_func.py
* This code uses the tabulated MESA EOS values to determine thermodynamic values of each particle instead of analytically solving for them
* This code has two functions
* The read_eos() function reads the entire sph.eos table and returns every read-in value
* The useeostable() uses the rho, u, and X (H1 composition) values to determine either the temperature, pressure, or mean molecular weight of each particle. For our purposes, we are finding the temperature and pressure of each particle

* It is important to note when to use each type of method. When running splot.py, you will be prompted to choose how you want the temperature/pressure values calculated. You should choose the analytic solution if you did not use the MESA EOS table in your simulation or if you do not have access to the MESA EOS table. You should use the MESA EOS table if you used that for your simulation. If possible, you should attempt to use the MESA EOS table for your simulation if you are planning on evolving the SPH remnant, as it accounts for other sources of energy including recombination/ionization and degeneracy pressures. This makes it more accurate to what MESA would expect when running the continued evolution

The following codes are used to determine which particles belong to the stellar remnant by creating a comp*.sph (for component) file.

## calccom.py

## component_best3.py

## component_reader.py
* This code reads in the new component file and saves the particle data for every particle that is bound, while omitting all other particles

## bound_particles.py
* This code uses the data from the component file along with the read in data and saves the particle data only for the bound particles to the star

--------------------------------------------------------------------------------------------------------------------------------------

# Expectations, Limitations, and Further MESA Analysis
1. The entropy file format might have an effect on the star's ability to converge down to a model. MESA's EOS tables work best when providing it a density and temperature (DT) pair. Otherwise, it will determine the temperature or density from another subroutine. If you are observing a star that is not in hydrostatic equilibrium, you may find that it is easiest to use a DT pair for your analysis
2. In your MESA inlist, you want to make sure that you include several variables in various settings:

### &star_job
    
    get_entropy_for_relax_from_eos = 'eosDE' !use eosDE for a density, specific thermal energy pair, eosDT for a density, temperature pair, or eosPT for a pressure, temperature pair
    timescale_for_relax_entropy = !this value should be smaller for stars in hydrostatic equilibrium (~1d-9), but should be increased for disrupted stars. If it is too high, it will skip portions of the star's evolution, and if it is too low, the star will run into errors converging
    
    relax_initial_composition = .true.
    relax_composition_filename = 'composition.dat'
    relax_initial_angular_momentum = .true.
    relax_angular_momentum_filename = 'angular_momentum.dat'
    relax_initial_entropy = .true.
    relax_entropy_filename = 'entropy.dat'

### &controls

    max_number_retries = 400 ! increase this value if your star fails after retrying too many times. This is especially important for stars not in hydrostatic equilibrium
    relax_max_number_retries = 500 ! increase this value if your star fails after retrying too many times in the entropy relaxation before the stellar evolution calculations begin
    
    initial_mass = ! this value is outputted when the files are generated and should be used for the resulting star. There will likely be a discrepancy between the MESA and StarSmasher value after several decimal places

3. As of right now, this code is not configured to allow for jump runs in the SPH, meaning all particles need to be present throughout the entire simulation. A good variable to help work around this would be the cc variable in the StarSmasher output, which saves parent-star data for the simulation
4. If you are having trouble using the eosDE method, you may have to edit a file in the MESA_DIR. The file can be found at

    $MESA_DIR/star/job/run_star_support.f90

and the line that needs to be edited is line 3062. It should read

    logT_guess = log10(min(T_guess_gas,T_guess_rad))

  5. As of now, this method has several limitations to the scenarios that it can apply towards.

## Symmetry

MESA assumes spherical symmetry, whereas StarSmasher and other SPH models can model 3D asymmetries. Because of this, it is important that your star from SPH is roughly spherical. The code should not have any issues creating the files, and MESA may be able to read them in. However, the resulting evolution will be that of a spherical star with spherical zones, even though the data used to define these zones are not spherical distributed and may not be indicative of the stellar structure (e.g. a tail or a "stretched out" star)

## Metallicity

  If you want to ensure that your SPH run is to be accurate to the MESA definitions of the thermodynamic values, you will need to use the MESA EOS tables. As of now, I only have access to tables that account for a metallicity of Z=0.02. This limitation can be simply (although tediously) remedied by writing an eos table file that accounts for a different metallicity values. However, this method should work if further MESA EOS tables were utilized assuming that they follow the desired format the the SPH uses to read the data in.
  
  This becomes increasingly complicated for older giant stars (i.e. Hydrogen-burning stars), whose metallicities may vary dramatically in the core of the star compared to its outer layers. If possible, the MESA EOS tables used would need to account for more than just one metallicity.
  
  It is important to note that the analytic EOS method is also a good estimate and may be an okay approximation for the thermodynamic values that would be used by MESA. As mentioned above, it is missing some factors in the pressure, energy, temperature, etc. that may offset the values and not provide as accurate of a stellar evolution as possible, but if the MESA EOS table is unable to work, the analytic EOS could be an alright way to estimate the evolution of the star for non-extreme models (where extreme could include things like very high-density models where degeneracy pressures would be important or very low density/temperature particles where recombination may significantly offset the energies and temperatures).


# References
* James Lombardi for assisting with the physical understanding of StarSmasher and its outputs, assisting with the writing of some of the codes, and for sharing prior Fortran codes that I was able to translate to Python for additional analysis
* Hans Vanderzyden for writing the initial Python files that read in the StarSmasher data and for writing the original splot.py code that I used as an interface
* Evan Bauer, Ebraheem Farag (on the mesa-users list) for assisting with using the files to relax the star and for helping alter the MESA code
