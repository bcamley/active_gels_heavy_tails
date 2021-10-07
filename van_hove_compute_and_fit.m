function [stable_fit,msd_out,xi_out,du] = van_hove_compute_and_fit(u,n,thin_by_nlag)

du = u - circshift(u,n); % Difference in u at two different lag times
du(1:n) = []; % First n values are useless...

msd_out = mean(du.^2);
xi_out = exp(mean(log(abs(du))));

try
    lastwarn('');
    if(thin_by_nlag)
        thin_by = floor(n*2);
        stable_fit = fitdist(du(1:thin_by:end)','stable');
    else
        stable_fit = fitdist(du(1:10:end)','stable');
    end
    [warnMsg, warnId] = lastwarn;
    betarange = paramci(stable_fit,'parameter','beta');
    beta_bad = max(abs(betarange))>0.5; % if we can't specify that beta = 0 better than this, fit was bad
    if(~isempty(warnMsg) || beta_bad)
        warning('Retrying with all data');
        stable_fit = fitdist(du(1:1:end)','stable');
    end
    
catch err
    getReport(err)
    stable_fit = NaN;
    alphas(i) = NaN;
    gams(i) = NaN;
end

