% Use this common function common to all models as a temporary substitute
% to reading the input file

function input_specs = InputRead_ABcat

input_specs.spec_names = {'A*','B*','*'};
%input_specs.N_0 = [30 10 60];                                                             % Initial state
input_specs.N_0 = [0 0 100];       

input_specs.stoich = [1 0 -1
    -1 0 1
    -1 1 0
    1 -1 0
    0 -1 1
    0 1 -1];         % stoichiometric matrix
input_specs.k = [4.72, 4.72/6.67, 2, 1, 0.4, 0];                                                         % Rate constants

% input_specs.stoich = [1 0 -1
%     -1 0 1
%     -1 1 0
%     1 -1 0
%     0 -1 1];         % stoichiometric matrix
% input_specs.k = [1, 1.5, 2, 1, 0.4];

input_specs.param_names = {};
for i = 1:length(input_specs.k)
    input_specs.param_names = [input_specs.param_names ['k_' num2str(i)]];
end

input_specs.t_final = 5;
input_specs.N_record = 1001;                                                                    % Arbitrary with no implications on the numerics, just used for visualizing average trajectories and checking it against the algebraic ODE solution
input_specs.fast_rxns = [1,2];
input_specs.eps = 1;                                                                        % Stiffness parameter

% Compute reactant matrix
input_specs.stoich_react = input_specs.stoich;
input_specs.stoich_react(input_specs.stoich_react > 0) = 0;
input_specs.stoich_react = -input_specs.stoich_react;

% Need to check that the choice of fast reactions is appropriate

input_specs.num_batches = 50;                   % Use 100 if you need even better sampling
input_specs.delta = 0.05;

end