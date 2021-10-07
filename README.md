# Code for the paper "Active gels, heavy tails, and the cytoskeleton"
### Written by Brian Camley and Daniel Swartz

To generate many of the plots and illustrate how the analysis works, the code "generate_default_plots.m" loads a file generated by the simulation and computes the van Hove distributions, fits them to the Levy form, and generates our plots. The name of the file analyzed is specified in generate_default_plots.m, and by default it analyzes "default_parameter_finite_dipole_long_0.mat" - this analysis will take a few minutes.

These files can be generated by run_default.m, which is called with a parameter simrun (just an integer to specify one of multiple simulation runs if necessary), dset (sets the dipole distance d), and "inputfilename" which sets the format of the output file as "filename_simrun" -- these can take a long time (12-24 hours depending on the computer in our experience; estimated times output). 

The main simulation code is in "tensor_sum_poisson_arbitrary_tensor.m" which is called with a structure of parameters.

To generate MSD exponents as a function of varying rheology (Fig. 4) use the command
[powerlaw_exps,betas,powerlaw_theory] = beta_vary_load_msd_compute('default_parameter_finite_dipole_long_0.mat');

To generate the Levy alphas when varying rheology (Fig. 10), use "alpha_as_function_of_beta_many.m" - this will require default_parameter_finite_dipole_long0-8.mat; these are large files and will be stored elsewhere.



Quasi-2D results use the Struve functions written by T.P. Theodoulidis (Theo2); these are included. See license_for_struve_functions.txt

