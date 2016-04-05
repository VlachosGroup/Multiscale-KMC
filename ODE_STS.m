% Orginary differential equation (ODE)
% Single time scale (STS)

% Computation and analysis are in the same script because it runs fast

function ODE_STS

clc; clear; fclose('all');

%% System specifications - read this from input file
spec_names = {'A*','B*','*'};
N_0 = [30, 10, 60];                  % initial condition 
t_final = 2;
epsilon = 10^-(3);
k = [1, 1.5, 2, 1, 0.4, 0];
k(1:2) = k(1:2) / epsilon;
S = [1 -1 -1 1 0 0
    0 0 1 -1 -1 1
    -1 1 0 0 1 -1]';            % rows: rxns, cols: species
[n_rxns, n_specs] = size(S);

S_react = [0 1 1 0 0 0
     0 0 0 1 1 0
     1 0 0 0 0 1]';

%% Solve the ODE system
[T,Y] = ode15s(@rate_eqns_param, [0, t_final], [N_0, zeros(1, n_specs*n_rxns)]);

    function dNdt = rate_eqns_param(t,z) 

        dNdt = zeros(n_specs + n_specs * n_rxns,1);     % account for species and sensitivities
        N = z(1:n_specs)';               % species populations
        C = z(n_specs+1:end);           % sensitivities
        
        [r, dr_dtheta, dr_dN] = rxn_rates(S_react, N, k);               % Get rate law information
        dNdt(1:n_specs) = S' * r;                                                           % Get differential changes in species populations
        
        % Indexes the matrix elements of dNdk into a vector for
        % differential equation
        % Try to vectorize this
        for i = 1:n_rxns                   
            dNdt(i*n_specs+1:(i+1)*n_specs) = S' * (dr_dN * C((i-1)*n_specs+1:i*n_specs) + dr_dtheta(:,i));
        end

    end

% Scale the sensitivities by epsilon, this saves us from having to code in
% additional functional dependence on epsilon
Y(:,4:9) = Y(:,4:9) / epsilon;
spec_pop_traj = Y(:,1:n_specs);
sens_traj = Y(:,n_specs+1:end);

%% Plot Results

% Plot species populations
figure
hold on
for spec = 1:n_specs
    plot(T,spec_pop_traj(:,spec))
end
hold off
box('on')
ax = gca;
ax.FontSize = 18;
title('Species Populations')
xlabel('time (s)','FontSize',18)
ylabel('population','FontSize',18)
legend('A','B','*');

% Plot sensitivities for each species
for spec = 1:n_specs
    
    inds = linspace(spec, n_rxns * n_specs - n_specs + spec, n_rxns);
    
    figure
    hold on
    plot(T, sens_traj(:,inds) )
    hold off
    box('on')
    xlabel('time (s)','FontSize',18)
    ylabel('sensitivity','FontSize',18)
    title([spec_names{spec} ' sensitivities'])
    ax = gca;
    ax.FontSize = 18;
    legend('k1','k2','k3','k4','k5','k6');
end


end