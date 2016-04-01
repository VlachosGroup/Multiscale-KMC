% This script simulates the following chemical reaction network

% 1. A(g) -> A*               k1
% 2.   A* -> A(g)             k2
% 3.   A* -> B*               k3
% 4.   B* -> A*               k4
% 5.   B* -> B(g)             k5

% All reactions are on the same time scale. Use this to analyze the
% effectiveness of the batch stopping criteria.

% The multi-scale sensitivity analysis of (M. Núñez and D.G. Vlachos, J. Chem. Phys. 142 (4), 044108 (2015)) is used.

clc; clear;
t_cpu_start = cputime;

% Change random seed to folder name
% currentDirectory = pwd;
% [upperPath, deepestFolder, ~] = fileparts(currentDirectory);
% rng(str2num(deepestFolder));

%% System Definition

% Inputs: Consider having all system specifications (parameters, stoichiometric matrix, N_initial, time scale? (can estimate from STS sim)) be read from an external file        
k = [1, 1.5, 2, 1, 0.4, 0];                                                         % Rate constants
stoich = [1 0 -1
    -1 0 1
    -1 1 0
    1 -1 0
    0 -1 1];                                                                        % stoichiometric matrix
N_initial = [30 60 10];                                                             % Initial state
%t_final = 4;
N_record = 1001;                                                                    % Arbitrary with no implications on the numerics, just used for visualizing average trajectories and checking it against the algebraic ODE solution 
fast_rxns = [1,2];
t_final = 100;
num_batches = 50;                   % Use 100 if you need even better sampling
delta = 0.05;

% System info
[n_params, n_specs] = size(stoich);
slow_rxns = linspace(1,n_params,n_params);
slow_rxns(fast_rxns) = [];                                                          % All the reactions which are not fast are slow
n_rxns_slow = length(slow_rxns);
n_rxns_fast = length(fast_rxns);
stoich_fast = stoich(fast_rxns,:);
stoich_slow = stoich(slow_rxns,:);
N_r = zeros(N_record,n_specs); 
N_int = zeros(1,n_specs);
N_int_r = zeros(N_record,n_specs);
N_int_prev = N_int;
t_r = linspace(0, t_final, N_record);                                               % Recording times, remains constant

%% Stochastic Simulation Loop

% Simulation initialization
N = N_initial;
N_prev = N;                                                                         % An initialization, will be truncated in taking statistics
t = 0;                                                                              % Initial macro time
N_int = zeros(1,n_specs);
ind_rec = 1;                                                                        % Keep track of which time point must be sampled
t_prev = 0;
n_events = 0;

% Sensitivity analysis Parameters
da_betabar_dtheta = zeros(n_rxns_slow,n_params);
W = zeros(1,n_params);
W_r = zeros(N_record,n_params);
W_prev = W;
dNbar_dtheta = zeros(n_params, n_specs);
dNbar_dtheta_r = zeros(n_params, n_specs, N_record);

while t < t_final                                                                   % (macro) Termination time controls the sampling
    
    disp(['Event # ' num2str(n_events)])
    disp(['Time: ' num2str(t) ' s'])
    
    % Record the current state as long as time >= t_sample
    while t >= t_r(ind_rec)
        
        % Record the species numbers
        N_r(ind_rec,:) = N_prev;
        N_int_r(ind_rec,:) = N_int_prev + N_prev * (t_r(ind_rec) - t_prev);
        W_r(ind_rec,:) = W_prev - sum(da_betabar_dtheta) * (t_r(ind_rec) - t_prev);
        dNbar_dtheta_r(:,:,ind_rec) = dNbar_dtheta;
        
        ind_rec = ind_rec + 1;                                                                              % Increment the recording index
    end
    
    % Visit micro scale to get fast-class averages and choice of macro
    % reaction
    batch_length_orig = 2;                                                          % Guess of the time-scale of the fast class
    [N_int_avg, N_origin, betabar, slow_rxn_to_fire, dt_macro, da_betabar_dtheta, dNbar_dtheta] = micro_scale(N, k, fast_rxns, slow_rxns, stoich_fast, batch_length_orig, num_batches, delta);
    N_prev = N_int_avg;
    N = N_origin + stoich_slow(slow_rxn_to_fire,:);                                 % which reaction will fire? And from which state?
    
    % Update time
    t_prev = t;
    t = t + dt_macro;
    n_events = n_events + 1;
    
    % Record previous values, must be done before you change them for the
    % next time step
    W_prev = W;
    dW = 1 / betabar(slow_rxn_to_fire) * da_betabar_dtheta(slow_rxn_to_fire,:) - sum(da_betabar_dtheta) * dt_macro;
    W = W + dW;
    
    % Integral values
    N_int_prev = N_int;
    N_int = N_int_prev + dt_macro * N_prev;
    
end

% Fill in the recording times that were missed
while ind_rec < N_record + 1
    
    N_r(ind_rec,:) = N_prev;
    N_int_r(ind_rec,:) = N_int_prev + N_prev * (t_r(ind_rec) - t_prev);
    W_r(ind_rec,:) = W_prev - sum(da_betabar_dtheta) * (t_r(ind_rec) - t_prev);
    dNbar_dtheta_r(:,:,ind_rec) = dNbar_dtheta;
    
    ind_rec = ind_rec + 1;                                          % Increment the recording index
end


disp('CPU time')
elapsed_cpu = cputime-t_cpu_start

%% Print Data into Output File

fidout = fopen('MSA_output.bin','w');
output_mat = [t_r', N_r, N_int_r, W_r];
fwrite(fidout,output_mat,'double');
fclose(fidout);

fidout = fopen('micro_derivs.bin','w');
fwrite(fidout,dNbar_dtheta_r,'double');
fclose(fidout);


