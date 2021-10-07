function [us,stablefit] = test_fit_2d(L,rho,Nits)

us = NaN*ones(1,Nits);
for s = 1:Nits
start_bound = poissrnd(round(rho*L^2));
xs = L*(rand(1,start_bound)-0.5);
ys = L*(rand(1,start_bound)-0.5);
rs = sqrt(xs.^2+ys.^2);

amplitudes = (rand(size(rs))-0.5);

umany = amplitudes./(rs);
us(s) = sum(umany);
end

stablefit = fitdist(us.','stable');