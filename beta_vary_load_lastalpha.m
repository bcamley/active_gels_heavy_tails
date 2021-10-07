% vary density
function [alphas,betas] = beta_vary_load_lastalpha(savefilename)
thin_by_nlag = false;
betas = linspace(0.01,0.98,15);
viscoelastic_convolve = true;
fontsize = 20;
measured_alpha_time = 100; % measure alpha at 100s lag time
%%
alphas = zeros(size(betas));
load(savefilename,'Hxtot','params','ts');
for j = 1:length(betas)
    params.beta = betas(j);
    if viscoelastic_convolve==true 
        u = convolve_fft_purepowerlaw(ts,Hxtot,params);
    else
        u = Hxtot / params.Gelastic;             % Treat as an elastic material with stiffness G(omega = 0)
    end
    u = real(u);
    n = round(measured_alpha_time/params.dt); % find the van hove, lag time = 100 s 
    [stable_fit,msd_out,xi_out,du] = van_hove_compute_and_fit(u,n,thin_by_nlag);
    %save(savefilename);
    alphas(j) = stable_fit.alpha;
    fprintf('beta = %3.3g, %d/%d \n',betas(j),j,length(betas));
end

% 
%% Plot the last-lag-time alpha as a function of beta
%figure

plot(betas,alphas,'o-','LineWidth',3);
%errorbar(betas,alphas_last,alphas_last-alphas_last_lower,alphas_last_upper-alphas_last,'o-','LineWidth',3);
set(gca,'FontSize',fontsize);
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
xlabel('Rheology exponent \beta')
ylabel(sprintf('van Hove \\alpha at t = %3.0f s',measured_alpha_time))