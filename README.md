# sph-to-mesa
This repository provides the code necessary to take an SPH model from StarSmasher and generate the necessary files to continue its evolution in MESA.

## Summary
Here is a description of the process and the following codes:
1. Ensure that the data from the original data mesa profile (i.e. profile15.data) contains all composition data for H1, He3, He4, C12, N14, O16, Ne20, and Mg24
2. Relax your star using StarSmasher. See https://jalombar.github.io/starsmasher/ for more information
3. In your working directory, run splot.py and run option 2. This will create a file called composition.sph
4. Run the SPH simulation and make sure that you know which output file(s) (out*.sph) that you will be analyzing
5. Run splot.py and use option 1. This will write 3 files to the directory: entropy.dat, composition.dat, and angular_momentum.dat
6. These three files can be used in the MESA test_suite relax_composition_j_entropy to continue the star's evolution

Here is a summary of each script:

## splot.py
### This code, originally created by Hans Vanderzyden and later modified by myself, acts as the main interface for the data analysis
### It's location should be in your working directory and it will import the other analysis files from their file locations
### There are currently two options for the code. Additional codes can be written for further analysis needs

The following scripts are for generating the necessary MESA files:

## readit.py
* This code, originally written by Hans Vanderzyden and later modified by myself, reads in the binary data from the out*.sph file(s) and saves all the data to a dictionary to be exported

## compsph_reaadit.py
* This code reads in composition data from the composition.sph file generated from the relaxed model
* The data is then saved as a list of lists which can be later read in in bestfit_total.py

## bestfit_total.py
* This code, translated from a Fortran routine by Dr. Jamie Lombardi, uses the data from readit.py and bins all the data based on their densities
* Following the binning process, the script determines various quantities including mass (m), radius (r), pressure (P), specific internal energy (u), density (rho), specific angular momentum (j), and composition of each bin
* These values are all written to a new file called bestfit.sph

## composition_j_entropy_reader.py
* This code reads in bestfit.sph from bestfit_total.py and creates a spline of the data, specifically of u vs x, rho vs x, j vs x, and composition vs x 
* Then it reads in a dummy file, angular_momentum.dat and assigns reasonable mass fraction (q) values to each spline
* Finally, it writes the data points of each spline and the corresponding q value to three new files, entropy.dat, composition.dat, and angular_momentum.dat
* The file is formatted to have the first line have the number of zones in the new star (and the number of elements in composition.dat) followed by a new line for every zone, starting with the q value, followed by the appropriate value to be analyzed. entropy.dat has density and specific energy (q   rho   u), angular_momentum.dat has the specific angular momentum (q   j), and composition.dat has the element abundances (q   H1   He3   He4   C12   N14   O16   Ne20   Mg24)

The following scripts are for generating the correct composition files:

## read_orig_comp.py
* This code reads in composition.dat, a file generated from an external script called composition.py, and saves the composition data for each zone (shell) of the star
* The code returns two numpy arrays. The first is the array of all q values that correspond with the shell, and the second is an 8D array of all 8 composition fractions of each shell.

## composition_spline.py
* This code uses the data from read_orig_comp.py and creates a spline of the data for each element
* This allows data to be interpolated between each data point corresponding to the original zones

## composition_fit.py
* This code uses composition_spline.py and the data from readit.py to now find the composition of each particle of the SPH model using a similar binning approach to what bestfit_total.py does
* It then writes all the data to composition.sph. Every line contains new data for each particle
