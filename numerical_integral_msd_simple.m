function msd = numerical_integral_msd_simple(beta,t)

msd = NaN*ones(size(t));

for j = 1:length(t)
    
    integrand = @(u) (1-cos(u)).*(abs(u).^(-2*beta))./(u.^2+t(j)^2);
    mic = 1e3; % max interval count
    abstol = 1e-4;
    reltol = 1e-3;
    msd(j) = (t(j)^(1+2*beta))*(quadgk(integrand,-Inf,-t(j),'MaxIntervalCount',mic,'AbsTol',abstol,'RelTol',reltol)+quadgk(integrand,-t(j),0,'MaxIntervalCount',mic,'AbsTol',abstol,'RelTol',reltol)+quadgk(integrand,0,t(j),'MaxIntervalCount',mic,'AbsTol',abstol,'RelTol',reltol)+quadgk(integrand,t(j),Inf,'MaxIntervalCount',mic,'AbsTol',abstol,'RelTol',reltol));    
end



end