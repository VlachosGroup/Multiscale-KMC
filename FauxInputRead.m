% Use this common function common to all models as a temporary substitute
% to reading the input file

function [spec_names, N_0, stoich, S_react, k, t_final, N_record, fast_rxns, eps] = FauxInputRead

spec_names = {'A*','B*','*'};
N_0 = [30 10 60];                                                             % Initial state

stoich = [1 0 -1
    -1 0 1
    -1 1 0
    1 -1 0
    0 -1 1
    0 1 -1];         % stoichiometric matrix
k = [1, 1.5, 2, 1, 0.4, 0];                                                         % Rate constants%

t_final = 4;
N_record = 1001;                                                                    % Arbitrary with no implications on the numerics, just used for visualizing average trajectories and checking it against the algebraic ODE solution
fast_rxns = [1,2];
eps = 0.01;                                                                        % Stiffness parameter

% Compute reactant matrix
S_react = stoich;
S_react(S_react > 0) = 0;
S_react = -S_react;

% Need to check that the choice of fast reactions is appropriate

end