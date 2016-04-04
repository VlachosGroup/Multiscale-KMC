function [rates,dr_dN,dr_dtheta] = rate_eqns(N, k)

% Reaction rates
rates = [k(1) * N(3)
    k(2) * N(1)
    k(3) * N(1)
    k(4) * N(2)
    k(5) * N(2)];

% dr_dN
dr_dN = [0 0 k(1)
    k(2) 0 0
    k(3) 0 0
    0 k(4) 0
    0 k(5) 0 ];

% dr_dtheta
dr_dtheta = diag([N(3) N(1) N(1) N(2) N(2)]);

end

