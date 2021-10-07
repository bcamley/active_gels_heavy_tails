% multiple toy stable fits
rng('shuffle');
RNGState = rng;

Ls = [5 10 20 40 80 160 400 800 1600];
rho = 0.1;
Nits = 5e4;

stablefits = cell(1,length(Ls));
alphas = NaN*ones(1,length(Ls));
alphas_upper = NaN*ones(1,length(Ls));
alphas_lower = NaN*ones(1,length(Ls));

for s = 1:length(Ls)
    [us,stablefit] = test_fit_2d(Ls(s),rho,Nits);
    stablefits{s} = stablefit;
    alphas(s) = stablefit.alpha;
    alpharange = paramci(stablefit,'parameter','alpha');
    alphas_lower(s) = alpharange(1); alphas_upper(s) = alpharange(2);
    fprintf('%d/%d: L = %3.1f alpha = %3.3f (CI: %3.3f, %3.3f)\n',s,length(Ls),Ls(s),alphas(s),alpharange(1),alpharange(2));
    clf
    subplot(1,2,1)
    errorbar(Ls,alphas,alphas-alphas_lower,alphas_upper-alphas,'o-','LineWidth',4);
    hold on
    plot(Ls,2*ones(size(Ls)),'--','LineWidth',3);
    subplot(1,2,2)
    errorbar(1./log(Ls),2-alphas,alphas_upper-alphas,alphas-alphas_lower,'s','LineWidth',4);
    hold on
    g = ~isnan(alphas);
    p = polyfit(1./log(Ls(g)),2-alphas(g),1);
    Linvf = linspace(0,max(1./log(Ls)),100);
    plot(Linvf,polyval(p,Linvf),'--','LineWidth',3);
    drawnow
end


% 
%% Final plot
fontsize = 20;

clf
subplot(1,2,1)

errorbar(Ls,alphas,alphas-alphas_lower,alphas_upper-alphas,'o-','LineWidth',4);
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
hold on
plot(Ls,2*ones(size(Ls)),'k--','LineWidth',3);
set(gca,'xscale','log');
xlabel('L'); ylabel('\alpha');
box off
subplot(1,2,2)
errorbar(1./log(Ls),2-alphas,alphas_upper-alphas,alphas-alphas_lower,'s','LineWidth',4);
hold on
g = ~isnan(alphas);
p = polyfit(1./log(Ls(g)),2-alphas(g),1);
Linvf = linspace(0,max(1./log(Ls)),100);
plot(Linvf,polyval(p,Linvf),'--','LineWidth',3);
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
xlabel('1/ln L'); ylabel('2-\alpha');
box off
%ylim([1 3]);