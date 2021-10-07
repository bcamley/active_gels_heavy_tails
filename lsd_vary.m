% vary Saffman-Delbruck length and density
rng('shuffle');
RNGState = rng;
set_default_parameters;
G2Ds = [10 100 1000];
%G2Ds = [0.01 0.1 1];
%G2Ds = 1e3;
krhos = [0.1 1 5]/params.tau;
params.tmax = params.tmax;
params.L = 80;
params.dim = 2;
params.dipole_resp_func = 'saffman_dipole_response';
savefilebasename = 'saffman_sweepG2D_rho_default';
viscoelastic_convolve = false;
fontsize = 20;
%params.cutoff = params.d;  %to create an explicit cutoff
j = 1;
for ii = 1:length(G2Ds)
    for kk = 1:length(krhos)
        params.G2D = G2Ds(ii);
        params.G0 = G2Ds(ii);
        params.Gelastic = params.G2D;
        params.kon_rhobulk = krhos(kk);
        savefilename = sprintf('%s%d.mat',savefilebasename,j);
        [Hxtot,ts,Numbound] = tensor_sum_poisson_arbitrary_tensor(params); save(savefilename);
        [u,tlag,xis,alphas,gams,msdfit,bincenters,binvalues,measured_stable_fit,dudump,fitdump,measured_normal_fit] = analysis_with_gamma(savefilename,'dontprintthis',viscoelastic_convolve,false);
        tlags{ii,kk} = tlag;
        alphas_all{ii,kk} = alphas;
        gams_all{ii,kk} = gams;
        fits_all{ii,kk} = measured_stable_fit;
        bins_all{ii,kk} = bincenters;
        values_all{ii,kk} = binvalues;
        normal_fits_all{ii,kk} = measured_normal_fit;
        close all
        fprintf('G2D = %3.3g, krho = %3.3g %d/%d \n',G2Ds(ii),krhos(kk),j,length(G2Ds)*length(krhos));
        j = j + 1;
    end
end

%%
%close all
subplot(length(G2Ds),length(krhos),1)
j = 1;
for kk = 1:length(krhos)
    for ii = 1:length(G2Ds)
        subplot(length(G2Ds),length(krhos),j)
        j = j + 1;
        fontsize = 14;
       
        plot(bins_all{ii,kk},values_all{ii,kk},'.','MarkerSize',12);
        hold on
        edges = bins_all{ii,kk}; % not actually edges, it's centers, but fine for getting our ranges
        plot_range = logspace(log10(0.5*min(abs(edges(edges>0)))),log10(2*max(abs(edges))),1e3);
        halpha = plot(plot_range,pdf(fits_all{ii,kk},plot_range),'--r','LineWidth',4);
        set(gca,'xscale','log','yscale','log','fontsize',fontsize);
        axis tight
        ylim_nonormal = ylim;
        plot(plot_range,pdf(normal_fits_all{ii,kk},plot_range),'.-','LineWidth',1);
        ylim(ylim_nonormal);
        legalpha = sprintf('\\alpha = %2.3f',fits_all{ii,kk}.alpha);
        legend(halpha,legalpha,'Location','southwest');
        %title(sprintf('G_{2D} = %3.1f, \\rho = %3.2f \\mum^{-2}: \\alpha = %2.2f',G2Ds(ii),krhos(kk)*params.tau,fits_all{ii,kk}.alpha));
        set(gca,'LineWidth',2,'TickLength',[0.025 0.025]);
        set(gca,'XTick',[1e-5 1e-4 1e-3 1e-2 1e-1])
        set(gca,'YTick',[1e-3 1 1e3]);
        
    end
end

%%
figure(2)
clf
 hold on
 fontsize = 20;
 for ii = 1:length(G2Ds)
     for kk = 1:length(krhos)
        plot(tlags{ii,kk},alphas_all{ii,kk},'LineWidth',3,'DisplayName',sprintf('G_{2D} = %3.2f, \\rho_{eff} = %3.2f \\mu m^{-2}',G2Ds(ii),krhos(kk)*params.tau));
     end
 end
set(gca,'xscale','log','yscale','linear','fontsize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
xlabel('Lag Time (sec)')
ylabel('Levy Stability Parameter \alpha')
legend