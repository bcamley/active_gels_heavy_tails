% initialize test for Poisson Levy
params = struct;
params.dim = 3;               % run in 3D
params.tmax = 50000;          % simulation time (s)
params.dt = 0.025;            %
params.kon_rhobulk = 0.1;     % kon_rhobulk * tau = average density; 
                              % kon_rhobulk * tau ~ 0.5 / micron^3 (Guo et
                              % al. Cell 2014)
params.L = 80;                % system size, microns
params.d = 0.4;               % dipole size, microns
params.tau = 5;               % mean dipole on time, from Guo et al. Cell 2014
params.F = 10;                % Dipole force, in pN (Guo et al. Cell 2014)
params.cutoff = 0;            % cutoff length: if cutoff > 0,
                              % exclude force dipoles from a distance
                              % cutoff of the origin
params.dipole_resp_func = ...
    'finite_dipole_response'; % options: finite_dipole_response ideal_dipole_response, saffman_dipole_response
                              
params.Gelastic = 10;         % Shear modulus assumed for purely elastic calculations (in Pa)
params.G0 = 0;                % Zero frequency shear modulus for viscoelastic case
params.beta = 0.17;           % rheological exponent
params.omega0 = 10;           % frequency, in rads/s
params.Gscale = 38;           % Shear modulus at frequency omega0, in Pa (from Hoffman et al. PNAS 2006)
params.Gint = 10;             % interior shear modulus for Saffman case, in Pa