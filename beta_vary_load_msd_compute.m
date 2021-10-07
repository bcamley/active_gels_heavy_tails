% vary density
function [powerlaw_exps,betas,powerlaw_theory] = beta_vary_load_msd_compute(savefilename)
thin_by_nlag = false;
betas = linspace(0.01,0.98,15);
viscoelastic_convolve = true;
fontsize = 20;
measured_alpha_time = 100; % measure alpha at 100s lag time
powerlaw_exps = zeros(size(betas));
load(savefilename,'Hxtot','params','ts');
for j = 1:length(betas)
    params.beta = betas(j);
    if viscoelastic_convolve==true 
        u = convolve_fft_purepowerlaw(ts,Hxtot,params);
    else
        u = Hxtot / params.Gelastic;             % Treat as an elastic material with stiffness G(omega = 0)
    end
    u = real(u);

    
    dt = params.dt;
    %tau = params.tau;
    tlag = logspace(-2,2); 
    nlag = unique(round(tlag/dt)); % converting to indicies
    tlag = dt*nlag; % lag times we will measure at
    msd = NaN*ones(size(nlag));
    i = 0;
    for n = nlag
       i = i + 1;
       du = u - circshift(u,n); % Difference in u at two different lag times
       du(1:n) = []; % First n values are useless... 
       msd(i) = mean(du.^2); 
    end
    
    g = (tlag < (params.tau/5)) & (tlag>0);
    msdfit = polyfit(log(tlag(g)),log(msd(g)),1);
    powerlaw_exps(j) = msdfit(1);
    
    %n = round(measured_alpha_time/params.dt); % find the van hove, lag time = 100 s 
    %[stable_fit,msd_out,xi_out,du] = van_hove_compute_and_fit(u,n,thin_by_nlag);
    %save(savefilename);
    %[u,tlag,xis,alphas,gams,msdfit] = 
    %[u,tlag,xis,alphas,gams,msdfit,bincenters,binvalues,measured_stable_fit,dudump,fitdump,measured_normal_fit,alphas_upper,alphas_lower] = analysis_with_gamma(savefilename,'dontprintthis',viscoelastic_convolve,false,betas(j));
    %alphas(j) = stable_fit.alpha;
    fprintf('beta = %3.3g, %d/%d \n',betas(j),j,length(betas));
end

load(savefilename,'params');
close all
clf
hold on
fontsize = 20;
plot(betas,powerlaw_exps,'ko-','LineWidth',4);
powerlaw_theory = sweep_beta_simple_msd(betas,params.tau);
hold on
plot(betas,powerlaw_theory,'--','LineWidth',3);
plot(betas,1+2*betas,'.-','LineWidth',3);
ylim([1 2]);
set(gca,'FontSize',fontsize);
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
xlabel('Rheology exponent \beta')
ylabel('MSD exponent \Delta')
legend('Simulation','Toy model','1+2\beta')
ylim([1 3]);