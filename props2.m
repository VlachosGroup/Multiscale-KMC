% Kinetic information on reaction rates/propensities along with their
% parameter derivatives
% Used in all model types

function [a, der] = props2(S_react, N,k)

% Compute a by going through the stoichiometry matrix and populations, use
% lines from your "PDE" code

    a = [k(1) * N(3)
         k(2) * N(1)
         k(3) * N(1)
         k(4) * N(2)
         k(5) * N(2)];      % Based on rate equations

    der = diag(a ./ k);
    % der = grad(a), the gradient of the rate equations
    % rows: reaction, cols: parameter
       
end