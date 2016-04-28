% Use this common function common to all models as a temporary substitute
% to reading the input file

function [spec_names, N_0, stoich, S_react, k, param_names, t_final, N_record, fast_rxns, eps, num_batches, delta] = FauxInputRead

spec_names = {'A','B','C'};
N_0 = [100 0 0];                                                             % Initial state

stoich = [-1 1 0
    1 -1 0
    0 -1 1];         % stoichiometric matrix
k = [1, 1.5, 2];                                                         % Rate constants%
param_names = {};
for i = 1:length(k)
    param_names = [param_names ['k_' num2str(i)]];
end

t_final = 5;
N_record = 1001;                                                                    % Arbitrary with no implications on the numerics, just used for visualizing average trajectories and checking it against the algebraic ODE solution
fast_rxns = [1,2];
eps = 0.01;                                                                        % Stiffness parameter

% TTS KMC only
num_batches = 50;                   % Use 100 if you need even better sampling
delta = 0.05;

% Compute reactant matrix
S_react = stoich;
S_react(S_react > 0) = 0;
S_react = -S_react;

end