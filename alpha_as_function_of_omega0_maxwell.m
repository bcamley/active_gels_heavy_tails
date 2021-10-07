clf
hold on

savefilenamebase = 'default_parameter_finite_dipole_long';
alphas_default = {};
for s = 0:8
    savefilename = sprintf('%s_%d',savefilenamebase,s);
    try
        [alphas,omega0s] = omega0_vary_load_lastalpha(savefilename);
        alphas_default{s+1} = alphas;
    catch err
        fprintf('Load or analysis failed for s = %d \n',s);
    end
    drawnow
end

savefilenamebase_d0p001 = 'default_d0p001';
alphas_d0p001 = {};
for s = 0:8
    savefilename = sprintf('%s_%d',savefilenamebase_d0p001,s);
    try
        [alphas,omega0s] = omega0_vary_load_lastalpha(savefilename);
        alphas_d0p001{s+1} = alphas;
    catch err
        fprintf('Load or analysis failed for s = %d \n',s);
    end
    drawnow
end

%% Collect the error bars
alphas_default_coll = cell2mat(alphas_default.');
alphas_default_mean = mean(alphas_default_coll);
alphas_default_err = std(alphas_default_coll)/sqrt(size(alphas_default_coll,1)); % standard error
errorbar(omega0s,alphas_default_mean,alphas_default_err,'ks-','LineWidth',3);

alphas_d0p001_coll = cell2mat(alphas_d0p001.');
alphas_d0p001_mean = mean(alphas_d0p001_coll);
alphas_d0p001_err = std(alphas_d0p001_coll)/sqrt(size(alphas_d0p001_coll,1)); % standard error
errorbar(omega0s,alphas_d0p001_mean,alphas_d0p001_err,'bs-','LineWidth',3);
set(gca,'xscale','log');