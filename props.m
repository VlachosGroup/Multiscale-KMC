% Kinetic information on all reaction evenets
% Could be used in a STS information or fed to micro scale of a TTS
% simulation (and parsed into pre-specified fast/slow reactions)

function [a, der] = props(N,k)

    a = [k(1) * N(3)
         k(2) * N(1)
         k(3) * N(1)
         k(4) * N(2)
         k(5) * N(2)];      % Based on rate equations

    der = [N(3)     0       0       0       0
           0        N(1)    0       0       0
           0        0       N(1)    0       0
           0        0       0       N(2)    0
           0        0       0       0       N(2)];              % Could generalize this by making it a n_params x n_reactions matrix for each parameter/reaction pair. Right now it assumed parameter = rate constant.
    % der = grad(a), the gradient of the rate equations
    % rows: reaction, cols: parameter
       
end