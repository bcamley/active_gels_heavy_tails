function [Hxtot,ts,Numbound] = tensor_sum_poisson_arbitrary_tensor(params)
%% Initializing parameters and setting up simulation
tmax = params.tmax; % maximum time of simulation in seconds
dt = params.dt; % time step in seconds
L = params.L; % Box size, box spans -L/2 to L/2
tau = params.tau; % sets our timescale -- the random switching time
ts = 0:dt:tmax-dt;%linspace(0,tmax,round(tmax/dt));
kon_rhobulk = params.kon_rhobulk;
cutoff = params.cutoff;%cutoff = 0;
% x-projection of trajectory
Hxtot = zeros(size(ts));
Numbound = zeros(size(ts));
if(~isfield(params,'dim'))
    params.dim = 3;
    warning('Assuming d = 3 dimensions');
end

switch(params.dim)
    case 3
        start_bound = round(kon_rhobulk*tau*L^3);
        max_numbound = round(4*kon_rhobulk*tau*L^3);
    case 2
        start_bound = round(kon_rhobulk*tau*L^2);
        max_numbound = round(4*kon_rhobulk*tau*L^2);
end
fprintf('Expecting a typical number of %d bound, current max is %d \n',start_bound,max_numbound);

isbound = false(1,max_numbound);
Hxi = zeros(1,max_numbound);
offtimes = NaN*ones(1,max_numbound);

% assume we start with start_bound bound
durations = -(tau)*log(rand(1,start_bound));

switch(params.dim)
    case 3
        xs = L*(rand(1,start_bound)-0.5);
        ys = L*(rand(1,start_bound)-0.5);
        zs = L*(rand(1,start_bound)-0.5);
        rs = sqrt(xs.^2+ys.^2+zs.^2);
    case 2
        xs = L*(rand(1,start_bound)-0.5);
        ys = L*(rand(1,start_bound)-0.5);
        rs = sqrt(xs.^2+ys.^2);
        zs = zeros(size(xs));
end

% This doesnt allow dipoles to get too close to the origin
if(cutoff>0)
    while any(rs < cutoff)
        num = sum(rs<cutoff);
        xs(rs<cutoff) = L*(rand(1,num)-0.5);
        ys(rs<cutoff) = L*(rand(1,num)-0.5);
        if(params.dim==3)
            zs(rs<cutoff) = L*(rand(1,num)-0.5);
            rs = sqrt(xs.^2+ys.^2+zs.^2);
        else
            rs = sqrt(xs.^2+ys.^2);
        end
    end
end

phis = 2*pi*rand(1,start_bound);
switch(params.dim)
    case 3
        thetas = acos((2*rand(1,start_bound)-1));
        bx = sin(thetas).*cos(phis);
        by = sin(thetas).*sin(phis);
        bz = cos(thetas);
    case 2
        bx = cos(phis);
        by = sin(phis);
        bz = zeros(size(phis));
end

isbound(1:start_bound) = true;
offtimes(1:start_bound) = durations;

% compute the initial displacement function H given these angles
%Hxi(isbound) = (F*d./(8*pi*rs.^2)).*(3*(bx.*xs./rs + by.*ys./rs + bz.*zs./rs).^2-1).*(xs./rs);
Hxi(isbound) = feval(params.dipole_resp_func,xs,ys,zs,rs,bx,by,bz,params); % using feval here gives a bit of a performance hit, but 
                                                                           % and can be sped up by hard-coding the
                                                                           % specific response function, but for clarity + reproducibility, not doing this
                                                                           % 
Hxtot(1) = sum(Hxi(isbound));
Numbound(1) = sum(isbound);
tic;

switch(params.dim)
    case 3
        mean_on_per_step = kon_rhobulk*dt*L^3
    case 2
        mean_on_per_step = kon_rhobulk*dt*L^2
end

%% Simulation Lööp
for s = 2:length(ts)
    % find those points which unbind this time step
    turnoff = (ts(s) >= offtimes) & (isbound);
    Hxtot(s) = Hxtot(s-1) - sum(Hxi(turnoff)); % remove the effect of all the force dipoles that turn off
    isbound(turnoff) = false;
    Hxi(turnoff) = 0; % should be redundant, ideally
    
    % compute the number of points who turn on:
    num_turnon = poissrnd(mean_on_per_step);
    
    if(num_turnon > 0)
        % compute the effect from turning these points on:
        xs = L*(rand(1,num_turnon)-0.5);
        ys = L*(rand(1,num_turnon)-0.5);
        if(params.dim==3)
            zs = L*(rand(1,num_turnon)-0.5);
            rs = sqrt(xs.^2+ys.^2+zs.^2);
        else
            zs = zeros(size(xs));
            rs = sqrt(xs.^2+ys.^2);
        end
        
        % This doesnt allow dipoles to get too close to the origin
        if(cutoff>0)
            while any(rs < cutoff)
                num = sum(rs<cutoff);
                xs(rs<cutoff) = L*(rand(1,num)-0.5);
                ys(rs<cutoff) = L*(rand(1,num)-0.5);
                if(params.dim ==3)
                    zs(rs<cutoff) = L*(rand(1,num)-0.5);
                    rs = sqrt(xs.^2+ys.^2+zs.^2);
                else
                    rs = sqrt(xs.^2+ys.^2);
                end
            end
        end
        
        phis = 2*pi*rand(1,num_turnon);
        switch(params.dim)
            case 3
                thetas = acos((2*rand(1,num_turnon)-1));
                bx = sin(thetas).*cos(phis);
                by = sin(thetas).*sin(phis);
                bz = cos(thetas);
            case 2
                bx = cos(phis);
                by = sin(phis);
                bz = zeros(size(phis));
        end
        offtimes_temp = ts(s)-(tau)*log(rand(1,num_turnon));
        Hxi_temp = feval(params.dipole_resp_func,xs,ys,zs,rs,bx,by,bz,params);
        
        % now need to assign which indices in the vector are unbound
        bind_inds = find(~isbound,num_turnon);
        if(length(bind_inds) < num_turnon)
            fprintf('Doubling the array sizes; this should not happen often! \n');
            Hxi = [Hxi zeros(size(Hxi))];
            offtimes = [offtimes zeros(size(offtimes))];
            isbound = [isbound false(size(isbound))];
            bind_inds = find(~isbound,num_turnon);
            if(length(bind_inds) < num_turnon)
                error('Rate of binding likely is unreasonably large')
            end
        end
        Hxi(bind_inds) = Hxi_temp;
        offtimes(bind_inds) = offtimes_temp;
        isbound(bind_inds) = true;
        Hxtot(s) = Hxtot(s) + sum(Hxi_temp);
    end
    Numbound(s) = sum(isbound);
    
    
    % Track progress...
    if(rem(s,1e5)==0)
        elapsedtime = toc;
        remaining_time_estimate = elapsedtime*((length(ts)-s)/s);
        fprintf('s = %d, %3.3g percent done, %3.f min elapsed, estimated %3.2f hr remaining \n',s,100*s/length(ts),elapsedtime/(60),remaining_time_estimate/(60*60));        
    end
    
    
end