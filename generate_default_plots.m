% generate most plots
data_file_name = 'default_parameter_finite_dipole_long_0.mat';
fig_direc = 'temp'
viscoelastic_convolve = true;
doprint = true;
[u,tlag,xis,alphas,gams,msdfit,bincenters,binvalues,measured_stable_fit,dudump,fitdump,measured_normal_fit,alphas_upper,alphas_lower] = analysis_with_gamma(data_file_name,fig_direc,viscoelastic_convolve,doprint);