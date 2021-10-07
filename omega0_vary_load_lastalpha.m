% vary density
function [alphas,omega0s] = omega0_vary_load_lastalpha(savefilename)
thin_by_nlag = false;
omega0s = logspace(-2,2,15);
viscoelastic_convolve = true;
fontsize = 20;
measured_alpha_time = 100; % measure alpha at 100s lag time
%%
alphas = zeros(size(omega0s));
load(savefilename,'Hxtot','params','ts');
for j = 1:length(omega0s)
    params.omega0 = omega0s(j);
    if viscoelastic_convolve==true 
        u = convolve_fft_maxwell_model(ts,Hxtot,params);
    else
        u = Hxtot / params.Gelastic;             % Treat as an elastic material with stiffness G(omega = 0)
    end
    u = real(u);
    n = round(measured_alpha_time/params.dt); % find the van hove, lag time = 100 s 
    [stable_fit,msd_out,xi_out,du] = van_hove_compute_and_fit(u,n,thin_by_nlag);
    %save(savefilename);
    %[u,tlag,xis,alphas,gams,msdfit] = 
    %[u,tlag,xis,alphas,gams,msdfit,bincenters,binvalues,measured_stable_fit,dudump,fitdump,measured_normal_fit,alphas_upper,alphas_lower] = analysis_with_gamma(savefilename,'dontprintthis',viscoelastic_convolve,false,betas(j));
    alphas(j) = stable_fit.alpha;
    fprintf('omega0 = %3.3g, %d/%d \n',omega0s(j),j,length(omega0s));
end

% 
%% Plot the last-lag-time alpha as a function of beta
%figure

plot(omega0s,alphas,'o-','LineWidth',3);
%errorbar(betas,alphas_last,alphas_last-alphas_last_lower,alphas_last_upper-alphas_last,'o-','LineWidth',3);
set(gca,'FontSize',fontsize);
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
xlabel('Maxwell model frequency \omega_0')
ylabel(sprintf('van Hove \\alpha at t = %3.0f s',measured_alpha_time))