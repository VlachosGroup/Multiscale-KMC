% Orginary differential equation (ODE)
% Single time scale (STS)

% Computation and analysis are in the same script because it runs fast

function ODE_STS

clc; clear; fclose('all');

% System specifications - read this from input file
addpath('../Network')
input_specs = InputRead_ABcat;
input_specs.k(input_specs.fast_rxns) = input_specs.k(input_specs.fast_rxns) / input_specs.eps;
[n_rxns, n_specs] = size(input_specs.stoich);

%% Solve the ODE system
[T,Y] = ode15s(@rate_eqns_param, linspace(0,input_specs.t_final,input_specs.N_record), [input_specs.N_0, zeros(1, n_specs*n_rxns)]);

    function dNdt = rate_eqns_param(t,z) 

        dNdt = zeros(n_specs + n_specs * n_rxns,1);     % account for species and sensitivities
        N = z(1:n_specs)';               % species populations
        C = z(n_specs+1:end);           % sensitivities
        
        [r, dr_dtheta, dr_dN] = rxn_rates(input_specs.stoich_react, N, input_specs.k);               % Get rate law information
        dNdt(1:n_specs) = input_specs.stoich' * r;                                                           % Get differential changes in species populations
        
        % Indexes the matrix elements of dNdk into a vector for
        % differential equation
        % Try to vectorize this
        for i = 1:n_rxns                   
            dNdt(i*n_specs+1:(i+1)*n_specs) = input_specs.stoich' * (dr_dN * C((i-1)*n_specs+1:i*n_specs) + dr_dtheta(:,i));
        end

    end

%% Plot Results

% Separate population and sensitivity data
spec_pop_traj = Y(:,1:n_specs);
sens_traj = Y(:,n_specs+1:end);

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
xlabel('time (s)','FontSize',18)
ylabel('spec. pop.','FontSize',18)
legend(input_specs.spec_names);

% Plot sensitivities for each species
for spec = 1:n_specs
    
    inds = linspace(spec, n_rxns * n_specs - n_specs + spec, n_rxns);
    sens = sens_traj(:,inds);
    sens(:,input_specs.fast_rxns) = sens(:,input_specs.fast_rxns) / input_specs.eps;        % Scale the sensitivities to eliminate dependence on eps
    
    
    figure
    hold on
    plot(T, sens )
    hold off
    box('on')
    xlabel('time (s)','FontSize',18)
    ylabel([input_specs.spec_names{spec} ' sensitivities'],'FontSize',18)
    ax = gca;
    ax.FontSize = 18;
    legend(input_specs.param_names);
end

disp('final time')
T(end)

disp('final species pops')
spec_pop_traj(end,:)

disp('final sensitivities')
sens_traj(end,:)

end