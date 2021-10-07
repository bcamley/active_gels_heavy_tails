% vary dipole size
rng('shuffle');
RNGState = rng;

set_default_parameters
ds = [0.001 0.1 0.4];

viscoelastic_convolve = false;
fontsize = 20;

tlags = {};
alphas_all = {};
gams_all = {};
Fd = params.F*params.d;

for j = 1:length(ds)
    params.d = ds(j);
    params.F = Fd/params.d; % keep F*d constant
    savefilename = sprintf('dipolesize_vary_keepFd_default_params%d.mat',j);
    [Hxtot,ts,Numbound] = tensor_sum_poisson_arbitrary_tensor(params); save(savefilename);
    [u,tlag,xis,alphas,gams] = analysis_with_gamma(savefilename,'dontprintthis',viscoelastic_convolve,false);
    tlags{j} = tlag;
    alphas_all{j} = alphas;
    gams_all{j} = gams;
end

%%
close all
clf
fontsize = 20;
subplot(1,2,1)
hold on

for j = 1:length(ds)
   plot(tlags{j},alphas_all{j},'LineWidth',3);
   legs{j} = sprintf('d = %3.3f \\mum',ds(j));
end

plot(tlags{j},1.5*ones(size(tlags{j})),'--','LineWidth',3);
legs{j+1} = 'Ideal dipole theory';
set(gca,'xscale','log','yscale','linear','fontsize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
xlabel('Lag Time (sec)')
ylabel('Levy Stability Parameter \alpha')
legend(legs)

subplot(1,2,2)
hold on
fontsize = 20;
for j = 1:length(ds)
   plot(tlags{j},gams_all{j},'LineWidth',3);
   legs{j} = sprintf('d = %3.4f \\mum',ds(j));
end
legs{j+1} = 'Ideal dipole theory';
Aconst = 0.0866371;
tlag = tlags{1};
gamma_theory = Aconst*(params.kon_rhobulk^(2/3))*params.d*(params.F/params.Gelastic)*(2*params.tau*(1-exp(-tlag/params.tau))).^(2/3);
plot(tlag,gamma_theory,'--','LineWidth',3);

set(gca,'xscale','log','yscale','linear','fontsize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
xlabel('Lag Time (sec)')
ylabel('Levy Scale Parameter \gamma [\mum]')
legend(legs)