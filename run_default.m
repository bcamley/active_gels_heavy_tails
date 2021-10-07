% default params
function RNGState = run_default(simrun,dset,inputfilename)
if(nargin<1)
    simrun = 1;
else
   try
     simrun = str2double(simrun);
   catch err
     getReport(err);
   end
end

simrun
rng('shuffle');
RNGState = rng;
rng(RNGState.Seed+857*simrun);
RNGState = rng
set_default_parameters
params.tmax = 5*params.tmax
%params.tmax = 0.001*params.tmax

if(nargin>1)
  try
    dset = str2double(dset);
  catch err
     getReport(err);
  end
  fprintf('Setting params.d = %3.3g \n',dset) ;
  params.d = dset;
end


if(nargin<3)
  savefilename = sprintf('default_parameter_finite_dipole_long_%d',simrun);
else
  savefilename = sprintf('%s_%d',inputfilename,simrun)
end

viscoelastic_convolve = true;
[Hxtot,ts,Numbound] = tensor_sum_poisson_arbitrary_tensor(params);
doprint = true;

figure_directory = sprintf('figs_%s',savefilename);
save(savefilename);
%%
%[u,tlag,xis,alphas,gams] = analysis_with_gamma(savefilename,figure_directory,viscoelastic_convolve,true);

end