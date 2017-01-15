% Differential algebraic equation (DAE), two time-scale ODE
% Two time scale (TTS)

function ODE_TTS

clc; clear; fclose('all');

% System inputs
input_specs = network_input();
[n_specs,n_rxns] = size(input_specs.stoich');

% Split the stoichiometry matrix into fast and slow components
Ss = input_specs.stoich';
Ss(:,input_specs.fast_rxns) = 0;
Sf = input_specs.stoich' - Ss;

% Transformation matrix for the variables to a new system
Ts = null(Sf','r')';
Tf = null(Ts,'r')';
T = [Tf; Ts];

[mf, ~] = size(Tf);
ms = n_specs - mf;
Tinv = inv(T);
Tf = T(1:mf,:);
Ts = T(mf+1:end,:);
fse = 2*mf + ms;                        % number of fast-scale equations

% Initial conditions
y_0 = T * input_specs.N_0';
Cy_0 = zeros(fse*n_rxns,1);

%% Solve the ODE system

% Create the mass matrix
num_eqns = n_specs + fse * n_rxns;
mass = eye(num_eqns);
mass(1:mf,:) = zeros(mf,num_eqns);                 % Algebraic equations for y_f
% Algebraic equations for C_yf
for j = 1:n_rxns
    mass(n_specs+1 + (j-1) * fse : n_specs + (j-1) * fse + 2*mf, :) = zeros(2*mf,num_eqns);
end

% Call ODE solver
options = odeset('Mass',mass);
[T,Y] = ode15s(@rate_eqns_param_TTS,linspace(0, input_specs.t_final, input_specs.N_record),[y_0; Cy_0],options);

    % Function for ODE solver
    function dz_dt = rate_eqns_param_TTS(t,z)
        
        dz_dt = zeros(num_eqns,1);
        
        % Get rate info
        [rates,dr_dtheta,dr_dN] = rxn_rates(input_specs.stoich_react, (Tinv * z(1:n_specs))', input_specs.k);
        
        % Differentials of y
        dz_dt(1:mf) = Tf * Sf * rates;              % fast
        dz_dt(mf+1:n_specs) = Ts * Ss * rates;            % slow
        
        % Differentials of Cy
        for i = 1:n_rxns
            dgdy  = Tf * Sf * dr_dN * Tinv(:,1:mf);
            dz_dt((i-1)*fse+1+n_specs : (i-1)*fse+n_specs+mf) = Tf * Sf * (dr_dN * Tinv * [z((i-1)*fse+1+n_specs :(i-1)*fse+n_specs+mf) + z((i-1)*fse+n_specs+mf+1:(i-1)*fse+n_specs+2*mf); z((i-1)*fse+n_specs+2*mf+1:i*fse+n_specs)]  + dr_dtheta(:,i));  % fast-fast and fast-slow        
            dz_dt((i-1)*fse+n_specs+mf+1 : (i-1)*fse+n_specs+2*mf) = z((i-1)*fse+1+n_specs :(i-1)*fse+n_specs+mf) + inv(dgdy) * Tf * Sf * dr_dtheta(:,i);                                                                                                                   % fast-fast
            dz_dt((i-1)*fse+n_specs+2*mf+1 : i*fse+n_specs) = Ts * Ss * (dr_dN * Tinv * [z((i-1)*fse+1+n_specs :(i-1)*fse+n_specs+mf) + z((i-1)*fse+n_specs+mf+1:(i-1)*fse+n_specs+2*mf); z((i-1)*fse+n_specs+2*mf+1:i*fse+n_specs)] + dr_dtheta(:,i));                                                                                                       % slow
        end
        
        
    end

% Convert back to oringinal variables
Tinv_aug = blkdiag(Tinv(:,1:mf),Tinv);                                                  % Add additional inverse lines because of the extra variables we have
Tbig = Tinv;
for h = 1:length(input_specs.k)
    Tbig = blkdiag(Tbig,Tinv_aug);
end
Yorig = Y * Tbig';
Y_sens = Yorig(:,n_specs+1:end);

%% Plots

% Plot species profiles
input_specs.plot_species_profiles(T,Yorig);

% Plot sensitivities for each species
for spec = 1:n_specs
    
    
    fast_sens = linspace(spec, 2 * n_specs * n_rxns - 2 * n_specs + spec ,n_rxns);
    slow_sens = linspace(spec + n_specs, 2 * n_specs * n_rxns - n_specs + spec,n_rxns);
    sens = Y_sens(:,fast_sens) + Y_sens(:,slow_sens);
    
    figure
    hold on
    for param = 1:n_rxns
        plot(T, sens(:,param))
    end
    
    hold off
    box('on')
    xlabel('Time (s)','FontSize',18)
    ylabel([input_specs.spec_names{spec} ' sensitivities'],'FontSize',18)
    ax = gca;
    ax.FontSize = 18;
    legend(input_specs.param_names);
    legend('boxoff')
end


end

