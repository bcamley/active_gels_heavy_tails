rerun_simulations = true;
rng('shuffle');
RNGState = rng;

set_default_parameters;
krhos = logspace(-2,0,4);
%krhos = linspace(1,10,4);

viscoelastic_convolve = true;
fontsize = 20;

for j = 1:length(krhos)
    params.kon_rhobulk = krhos(j);
    savefilename = sprintf('density_vary_default%d.mat',j);
    if(rerun_simulations)
        [Hxtot,ts,Numbound] = tensor_sum_poisson_arbitrary_tensor(params); save(savefilename);
    end
    
    [u,tlag,xis,alphas,gams] = analysis_with_gamma(savefilename,'dontprintthis',viscoelastic_convolve,false);
    tlags{j} = tlag;
    alphas_all{j} = alphas;
    gams_all{j} = gams;
    close all
    fprintf('Density at kbulkrho = %3.3g, %d/%d \n',krhos(j),j,length(krhos));
end

%%
close all
clf
hold on
fontsize = 20;
for j = 1:length(krhos)
   plot(tlags{j},alphas_all{j},'LineWidth',3);
   legs{j} = sprintf('\\rho_{avg} = %3.2f \\mum^{-3}',krhos(j)*(params.tau));
end
set(gca,'xscale','log','yscale','linear','fontsize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
xlabel('Lag Time (sec)')
ylabel('Levy Stability Parameter \alpha')
legend(legs)