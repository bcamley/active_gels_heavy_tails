fontsize=20;
clf
t = linspace(0,50,500);
tstart = 15;
ton = 5;
F = (t>tstart) & (t<(tstart+ton));

clf
hold on
set(gca,'FontSize',fontsize)
set(gca,'LineWidth',2,'TickLength',[0.025 0.025])
uth = @(beta,t) (1/(pi))*(gamma(1-beta)*sin(pi*beta)/beta)*( ((t-tstart).^beta).*(t>tstart) - (t>(tstart+ton)).*(t-tstart-ton).^beta);
uth_nopref = @(beta,t) ( ((t-tstart).^beta).*(t>tstart) - (t>(tstart+ton)).*(t-tstart-ton).^beta);
%plot(t,F);
h1=plot(t,uth_nopref(0,t),'LineWidth',8);
h2=plot(t,uth(0.15,t),'LineWidth',12);
h3=plot(t,0.3*uth_nopref(1,t),'LineWidth',5);
xlim([10 30]);
xlabel('Time')
ylabel('u^{(1)}(t)')
%set(gca,'FontSize',48)
legend('\beta = 0','\beta = 0.15','\beta = 1');