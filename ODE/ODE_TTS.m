% Differential algebraic equation (DAE), two time-scale ODE
% Two time scale (TTS)

function ODE_TTS

clc; clear; fclose('all');

% System inputs
addpath('../Network')
input_specs = InputRead_ABcat;
[M,R] = size(input_specs.stoich');

% Split the stoichiometry matrix into fast and slow components
Ss = input_specs.stoich';
Ss(:,input_specs.fast_rxns) = 0;
Sf = input_specs.stoich' - Ss;

% Transformation matrix for the variables to a new system
Ts = null(Sf','r')';
Tf = null(Ts,'r')';
T = [Tf; Ts];

[mf, ~] = size(Tf);
ms = M - mf;
Tinv = inv(T);
Tf = T(1:mf,:);
Ts = T(mf+1:end,:);
fse = 2*mf + ms;                        % number of fast-scale equations

% Initial conditions
y_0 = T * input_specs.N_0';
Cy_0 = zeros(fse*R,1);

%% Solve the ODE system

% Create the mass matrix
num_eqns = M + fse * R;
mass = eye(num_eqns);
mass(1:mf,:) = zeros(mf,num_eqns);                 % Algebraic equations for y_f
% Algebraic equations for C_yf
for j = 1:R
    mass(M+1 + (j-1) * fse : M + (j-1) * fse + 2*mf, :) = zeros(2*mf,num_eqns);
end

% Call ODE solver
options = odeset('Mass',mass);
[T,Y] = ode15s(@rate_eqns_param_TTS,linspace(0, input_specs.t_final, input_specs.N_record),[y_0; Cy_0],options);

    % Function for ODE solver
    function dz_dt = rate_eqns_param_TTS(t,z)
        
        dz_dt = zeros(num_eqns,1);
        
        % Get rate info
        [rates,dr_dtheta,dr_dN] = rxn_rates(input_specs.stoich_react, (Tinv * z(1:M))', input_specs.k);
        
        % Differentials of y
        dz_dt(1:mf) = Tf * Sf * rates;              % fast
        dz_dt(mf+1:M) = Ts * Ss * rates;            % slow
        
        % Differentials of Cy
        for i = 1:R
            dgdy  = Tf * Sf * dr_dN * Tinv(:,1:mf);
            dz_dt((i-1)*fse+1+M : (i-1)*fse+M+mf) = Tf * Sf * (dr_dN * Tinv * [z((i-1)*fse+1+M :(i-1)*fse+M+mf) + z((i-1)*fse+M+mf+1:(i-1)*fse+M+2*mf); z((i-1)*fse+M+2*mf+1:i*fse+M)]  + dr_dtheta(:,i));  % fast-fast and fast-slow        
            dz_dt((i-1)*fse+M+mf+1 : (i-1)*fse+M+2*mf) = z((i-1)*fse+1+M :(i-1)*fse+M+mf) + inv(dgdy) * Tf * Sf * dr_dtheta(:,i);                                                                                                                   % fast-fast
            dz_dt((i-1)*fse+M+2*mf+1 : i*fse+M) = Ts * Ss * (dr_dN * Tinv * [z((i-1)*fse+1+M :(i-1)*fse+M+mf) + z((i-1)*fse+M+mf+1:(i-1)*fse+M+2*mf); z((i-1)*fse+M+2*mf+1:i*fse+M)] + dr_dtheta(:,i));                                                                                                       % slow
        end
        
        
    end

% Convert back to oringinal variables
Tinv_aug = blkdiag(Tinv(:,1:mf),Tinv);                                                  % Add additional inverse lines because of the extra variables we have
Tbig = Tinv;
for h = 1:length(input_specs.k)
    Tbig = blkdiag(Tbig,Tinv_aug);
end
Yorig = Y * Tbig';
Y_sens = Yorig(:,M+1:end);

%% Plots
set(0,'defaultlinelinewidth',1.5)
set(0,'defaultaxeslinewidth',2)
ticklabelfontsize = 20;
axislabelfontsize = 24;


% Plot species populations
figure
hold on
for spec = 1:M
    plot(T,Yorig(:,spec))
end
hold off
box('on')
ax = gca;
ax.FontSize = 18;
xlabel('time (s)','FontSize',18)
ylabel('spec. pop.','FontSize',18)
legend(input_specs.spec_names);

% Plot sensitivities for each species
for spec = 1:M
    
    
    fast_sens = linspace(spec, 2 * M * R - 2 * M + spec ,R);
    slow_sens = linspace(spec + M, 2 * M * R - M + spec,R);
    sens = Y_sens(:,fast_sens) + Y_sens(:,slow_sens);
    
    figure
    hold on
    for param = 1:R
        plot(T, sens(:,param))
    end
    
    hold off
    box('on')
    xlabel('time (s)','FontSize',18)
    ylabel([input_specs.spec_names{spec} ' sensitivities'],'FontSize',18)
    ax = gca;
    ax.FontSize = 18;
    legend(input_specs.param_names);
    legend('boxoff')
end

%% Output file


end

