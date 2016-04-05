% Kinetic information on reaction rates/propensities along with their
% parameter and species derivatives

function [rxn_rates, dr_dk, dr_dN] = rxn_rates(S_rcnt, spec_pops, rate_conts)

% Inputs
% S_rcnt: n_rxns x n_species
% spec_pops: 1 x n_species
% rate_conts: 1 x n_rxns

% Outputs
% rxn_rates: n_rxns x 1
% dr_dk: n_rxns x n_rxns 
% dr_dN: n_rxns x n_species

n_rxns = length(rate_conts);
n_species = length(spec_pops);

% Compute rate constants/propensities
rxn_rates = rate_conts' .* prod(repmat(spec_pops, n_rxns, 1) .^ S_rcnt, 2);           % Uses reaction rate formula: r = k * N1 ^ n1 * N2 ^ n2 * ...

% Compute gradient of reaction rates with respect to rate constants
dr_dk = diag(prod(repmat(spec_pops, n_rxns, 1) .^ S_rcnt, 2));

% Compute gradient of reaction rates with respect to species populations
% Should find a faster, vectorized way of doing this
% Also put a flag so it only evaluates in ODE mode
dr_dN = zeros(n_rxns,n_species);
for i = 1:n_rxns
    for j = 1:n_species
        if S_rcnt(i,j) == 0
            dr_dN(i,j) = 0;
        else
            vector = zeros(1,n_species);
            vector(j) = 1;
            dr_dN(i,j) = S_rcnt(i,j) * rate_conts(i) * prod(spec_pops .^ (S_rcnt(i,:) - vector) );
        end
    end
end

end