% plot exponent vs beta
function Deltas = sweep_beta_simple_msd(betas,tau)
%betas = linspace(0,1,25);
t = logspace(-3,3,500);
msds = {};
Deltas = NaN*ones(size(betas));
for k = 1:length(betas)
    msd = numerical_integral_msd_simple(betas(k),t);
    msds{k} = msd;
    p = polyfit(log(t(t<0.1/tau)),log(msd(t<0.1/tau)),1);
    Deltas(k) = p(1);
end