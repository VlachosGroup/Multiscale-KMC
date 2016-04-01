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
k = [1, 1.5, 2, 1, 0.4];                                                         % Rate constants
stoich = [1 0 -1
    -1 0 1
    -1 1 0
    1 -1 0
    0 -1 1];                                                                        % stoichiometric matrix
N_initial = [30 10 60];                                                             % Initial state
%t_final = 4;
N_record = 1001;                                                                    % Arbitrary with no implications on the numerics, just used for visualizing average trajectories and checking it against the algebraic ODE solution 
t_final = 100;
fast_rxns = [1,2];
eps = 0.001;                                                                        % Stiffness parameter
k(fast_rxns) = k(fast_rxns) / eps;

% System info
[n_params, n_specs] = size(stoich);
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
da_dtheta = zeros(n_params,n_params);
W = zeros(1,n_params);
W_r = zeros(N_record,n_params);
W_prev = W;

while t < t_final                                                                   % (macro) Termination time controls the sampling
    
    disp(['Event # ' num2str(n_events)])
    disp(['Time: ' num2str(t) ' s'])
    
    % Record the current state as long as time >= t_sample
    while t >= t_r(ind_rec)                                                     % If you record after you compute the next step, but before you update data, then you won't need all the prev variables
        
        % Record the species numbers
        N_r(ind_rec,:) = N_prev;
        N_int_r(ind_rec,:) = N_int_prev + N_prev * (t_r(ind_rec) - t_prev);
        W_r(ind_rec,:) = W_prev - sum(da_dtheta) * (t_r(ind_rec) - t_prev);
        
        ind_rec = ind_rec + 1;                                                                              % Increment the recording index
    end
        
    [a, da_dtheta] = props(N,k);
    a_sum = sum(a);
    rxn_to_fire = min(find(rand(1)<cumsum(a/sum(a_sum))));     
    N_prev = N;
    N = N + stoich(rxn_to_fire,:);                                 % which reaction will fire? And from which state?
    
    % Update time
    dt = log(1/rand(1))/a_sum;
    t_prev = t;
    t = t + dt;
    n_events = n_events + 1;
    
    % Record previous values, must be done before you change them for the
    % next time step
    W_prev = W;
    dW = 1 / a(rxn_to_fire) * da_dtheta(rxn_to_fire,:) - sum(da_dtheta) * dt;
    W = W + dW;
    
    % Integral values
    N_int_prev = N_int;
    N_int = N_int_prev + dt * N_prev;
    
end

% Fill in the recording times that were missed
while ind_rec < N_record + 1
    
    N_r(ind_rec,:) = N_prev;
    N_int_r(ind_rec,:) = N_int_prev + N_prev * (t_r(ind_rec) - t_prev);
    W_r(ind_rec,:) = W_prev - sum(da_dtheta) * (t_r(ind_rec) - t_prev);
    
    ind_rec = ind_rec + 1;                                          % Increment the recording index
end


disp('CPU time')
elapsed_cpu = cputime-t_cpu_start

%% Print Data into Output File

fidout = fopen('MSA_output.bin','w');
output_mat = [t_r', N_r, N_int_r, W_r];
fwrite(fidout,output_mat,'double');
fclose(fidout);