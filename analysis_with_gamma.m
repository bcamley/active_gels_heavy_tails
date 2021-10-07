%% Load in data and set up parameters
function [u,tlag,xis,alphas,gams,msdfit,bincenters,binvalues,measured_stable_fit,dudump,fitdump,measured_normal_fit,alphas_upper,alphas_lower] = analysis_with_gamma(data_file_name,fig_direc,viscoelastic_convolve,doprint,beta,truncu)
%clear



thin_by_nlag = false;

dudump = NaN;
fitdump = NaN;
if(nargin<4)
    doprint = false;
end
load(data_file_name,'Hxtot','params','ts')

if(~isfield(params,'truncu'))
    truncate = false;
else
    truncate = true;
    if(nargin<6)
        truncu = params.truncu;
    end
end

if(nargin<5)
    beta = params.beta;
else
    params.beta = beta;
end

%load('simulation_data.mat')
try
    mkdir(fig_direc)
catch
    fprintf('Writing to existing figure directory %s',fig_direc)
    pause(1)
end

%u = Hxtot;
%%
if viscoelastic_convolve==true 
    u = convolve_fft_purepowerlaw(ts,Hxtot,params);
    % warning('doing pure power law');
    %u = convolve_fft(ts,u,params); % Takes data from elastic -> viscoelastic case
else
    u = Hxtot / params.Gelastic;             % Treat as an elastic material with stiffness G(omega = 0)
end
u = real(u);

if(truncate)
   %u = truncu*tanh(u./truncu); % will max at +/- truncu
   u0 = u;
   u(u0>truncu) = truncu + (u0(u0>truncu)-truncu)/2;
   u(u0<-truncu) = -truncu + (u0(u0<-truncu)+truncu)/2;
end

if(isfield(params,'localization_noise'))
   u = u + (params.localization_noise)*randn(size(u)); 
end

figure(15)
plot(u)
    

dt = params.dt;
tau = params.tau;
tlag = logspace(-2,2); 
nlag = unique(round(tlag/dt)); % converting to indicies
tlag = dt*nlag; % lag times we will measure at
lag_time_van_hove = 2*tau; % lag time at which we compute the Van Hove Correlation Plot

msd = zeros(size(nlag)); % Mean squared displacement at each lag time
xis = zeros(size(nlag)); % Rescaling factor as a function of t ~ gamma
gams = zeros(size(nlag)); % Value of gamma at each lag time
alphas = zeros(size(nlag)); % Value of levu alpha at each lag time 
alphas_upper = zeros(size(nlag));
alphas_lower = zeros(size(nlag));

fontsize = 20;
%% Make all the Figures

% Plot for a sample trajectory
figure(1)
hold on
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
title('Sample Trajectory Projected Along x-Axis')
xlabel('Time (sec)')
ylabel('Displacement u_x (μm)')
plot(ts(1:10:10000),u(1:10:10000),'LineWidth',2)

% Plot of MSD
figure(2)
hold on
set(gca,'FontSize',fontsize,'xscale','log','yscale','log')
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
xlabel('Lag Time (sec)')
ylabel('MSD (μm^2)')
xline(tau,'--b','\tau','LineWidth',3,'FontSize',24);


% Van Hove plot
figure(3)
hold on
set(gca,'xscale','log','yscale','log','fontsize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
title(['Van Hove Correlation Plot: Lag Time = ' num2str(lag_time_van_hove,'%.2f') 's'])
xlabel('Displacement \Deltau_x (μm)')
ylabel('Probability Density P(\Deltau_x)')
legend

% Plot of xi vs time
figure(4)
hold on
set(gca,'xscale','log','yscale','log','fontsize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
xlabel('Lag Time (sec)')
ylabel('Rescaling Parameter \xi')

% Plot of rescaled probabilities
figure(5)
hold on
set(gca,'xscale','log','yscale','log','fontsize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
%title('Rescaled Van Hove Correlation Plots')
xlabel('Displacement \Deltau_x/\xi ')
ylabel('Probability Density P(\Deltau_x/\xi)')
legend

% Plot of alpha vs time
figure(6)
hold on
set(gca,'xscale','log','yscale','linear','fontsize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
xlabel('Lag Time (sec)')
ylabel('Lévy \alpha')

figure(7)
hold on
set(gca,'xscale','log','yscale','linear','fontsize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
xlabel('Lag Time (sec)')
ylabel('Levy Width Parameter \gamma')

%% Actaully Calculating the curves to plot
i = 0;
for n = nlag
    i = i+1;
    
    [stable_fit,msd_out,xi_out,du] = van_hove_compute_and_fit(u,n,thin_by_nlag);
    msd(i) = msd_out;
    xis(i) = xi_out;
    alphas(i) = stable_fit.alpha;
    alpha_range = paramci(stable_fit,'parameter','alpha');
    alphas_upper(i) = alpha_range(2);
    alphas_lower(i) = alpha_range(1);
    gams(i) = stable_fit.gam;
    if(alphas(i)>1.95)
        dudump = du;
        fitdump = stable_fit;
    end

       
    if mod(i,10) == 0
        figure(5)
        [centers,values,plot_range] = log_histogram_pdf(du/xis(i));
        plot(centers,values,'.','MarkerSize',18,'DisplayName',['lag time = ' num2str(dt*n,'%.2f') 's'])
        %histogram(du/xis(i),'Normalization','pdf','Displaystyle','stairs','lineWidth',3,'DisplayName',['lag time = ' num2str(dt*n) 's'])
        
        if(i == 40)
            xl = xlim;
            plot_range = logspace(log10(xl(1)),log10(xl(2)),1e3);
            try
                stable_fit_rescale = fitdist(du(1:10:end)'/xis(i),'stable');
                plot(plot_range,pdf(stable_fit_rescale,plot_range),'--r','LineWidth',4,'DisplayName',['Fit with \alpha = ' num2str(stable_fit_rescale.alpha,'%.2f')]);
            catch err
                getReport(err)
                stable_fit_rescale = NaN; 
            end
                
        end
    end
    
end

figure(1)
if(doprint)
if viscoelastic_convolve == false
    saveas(gca,sprintf('%s/sample_trajectory.pdf',fig_direc))
else
    saveas(gca,sprintf('%s/visco_sample_trajectory.pdf',fig_direc))
end
end

figure(2)
g = (tlag < (params.tau/5)) & (tlag>0);
msdfit = polyfit(log(tlag(g)),log(msd(g)),1);
pmsd = plot(tlag,msd,'ko-','LineWidth',4,'MarkerSize',10,'DisplayName','Simulation');
pfit = plot(tlag,exp(polyval(msdfit,log(tlag))),'Linewidth',4,'DisplayName',['Power Law With Exponent ' num2str(msdfit(1),'%.2f')]);
legend([pfit,pmsd],['Best fit to t^\Delta with \Delta = ' num2str(msdfit(1),'%.2f')],'Simulation');

if(doprint)
if viscoelastic_convolve == false
    saveas(gca,sprintf('%s/MSD.pdf',fig_direc))
else
    saveas(gca,sprintf('%s/viscoMSD.pdf',fig_direc))
end
end

figure(3)
hold on
n_van_hove = round(lag_time_van_hove/dt);
du = u - circshift(u,n_van_hove);
if(thin_by_nlag)
    thin_by = floor(n_van_hove*2)
    stabledist = fitdist(du(1:thin_by:end)','stable');
else
    stabledist = fitdist(du(1:10:end)','stable');
end

%stabledist = fitdist(du(1:10:end)','stable');
%[~,edges] = histcounts(log10(du(du>0)));
%xlim([plot_range(1) plot_range(end)])
[centers,values,plot_range] = log_histogram_pdf(du);
bincenters = centers;
binvalues = values;
measured_stable_fit = stabledist;
measured_normal_fit = fitdist(du(1:10:end)','normal');
plot(centers,values,'.','MarkerSize',18,'DisplayName','Trajectory Data');
%histogram(du,plot_range,'Displaystyle','stairs','lineWidth',3,'Normalization','pdf','DisplayName','Trajectory Data')
plot(plot_range,pdf(stabledist,plot_range),'--r','LineWidth',4,'DisplayName',['Fit with \alpha = ' num2str(stabledist.alpha,'%.2f')]);

if(doprint)
if viscoelastic_convolve == false
    saveas(gca,sprintf('%s/vanhovelag.pdf',fig_direc))
else
    saveas(gca,sprintf('%s/viscovanhovelag.pdf',fig_direc))
end
end

figure(4)
plot(tlag,xis,'LineWidth',3)
if(doprint)
if viscoelastic_convolve == false
    saveas(gca,sprintf('%s/xivstime.pdf',fig_direc))
else
    saveas(gca,sprintf('%s/viscoxivstime.pdf',fig_direc))
end
end

figure(5)
if(doprint)
if viscoelastic_convolve == false
    saveas(gca,sprintf('%s/rescaledprobs.pdf',fig_direc))
else
    saveas(gca,sprintf('%s/viscorescaledprobs.pdf',fig_direc))
end
end

figure(6)
%plot(tlag,alphas,'o-','LineWidth',4,'MarkerSize',10)
errorbar(tlag,alphas,alphas-alphas_lower,alphas_upper-alphas,'o-','LineWidth',4,'MarkerSize',10);
if(doprint)
if viscoelastic_convolve == false
    saveas(gca,sprintf('%s/alphavstime.pdf',fig_direc));
else
    saveas(gca,sprintf('%s/viscoalphavstime.pdf',fig_direc))
end
end

figure(7)
plot(tlag,gams,'s-','LineWidth',4,'MarkerSize',10)
Aconst = 0.0866371;

gamma_theory = Aconst*(params.kon_rhobulk^(2/3))*params.d*(params.F/params.Gelastic)*(2*params.tau*(1-exp(-tlag/params.tau))).^(2/3);
plot(tlag,gamma_theory,'--','LineWidth',3);
if(doprint)
if viscoelastic_convolve == false
    saveas(gca,sprintf('%s/gammavstime.pdf',fig_direc))
else
    saveas(gca,sprintf('%s/viscogammavstime.pdf',fig_direc))
end
end

