% The multi-scale sensitivity analysis of M. Núñez and D.G. Vlachos, J. Chem. Phys. 142 (4), 044108 (2015) is used.

clc; clear; fclose('all');  t_cpu_start = cputime;

%% User Input - change to read everything from the input files

rng(12345);     % Set random seed
input_specs = InputRead_ABcat;
out_file = fopen('MSA_debug.txt','w');

%% Stochastic Simulation Loop

% System Information
[n_params, n_specs] = size(input_specs.stoich);
slow_rxns = linspace(1, n_params, n_params);
slow_rxns(input_specs.fast_rxns) = [];                                                          % All the reactions which are not fast are slow
n_rxns_slow = length(slow_rxns);
n_rxns_fast = length(input_specs.fast_rxns);
stoich_fast = input_specs.stoich(input_specs.fast_rxns,:);
stoich_slow = input_specs.stoich(slow_rxns,:);
N_r = zeros(input_specs.N_record,n_specs); 
N_int = zeros(1,n_specs);
N_int_r = zeros(input_specs.N_record,n_specs);
N_int_prev = N_int;
t_r = linspace(0, input_specs.t_final, input_specs.N_record);                                               % Recording times, remains constant

% Simulation initialization
N = input_specs.N_0;
N_prev = N;                                                                         % An initialization, will be truncated in taking statistics
t = 0;                                                                              % Initial macro time
N_int = zeros(1,n_specs);
ind_rec = 1;                                                                        % Keep track of which time point must be sampled
t_prev = 0;
n_events = 0;

% Sensitivity analysis Parameters
da_betabar_dtheta = zeros(n_rxns_slow,n_params);
W = zeros(1,n_params);
W_r = zeros(input_specs.N_record,n_params);
W_prev = W;
dNbar_dtheta = zeros(n_params, n_specs);
dNbar_dtheta_r = zeros(n_params, n_specs, input_specs.N_record);

while t < input_specs.t_final                                                                   % (macro) Termination time controls the sampling
    
    
    
    % Record the current state as long as time >= t_sample
    while t >= t_r(ind_rec)
        
        fprintf(out_file, 'Macroscopic steps executed: %d\n', n_events);               % Eventually write this into an output file
        fprintf(out_file, 'KMC time: %d s\n\n', t);
        
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
    [N_int_avg, N_origin, betabar, slow_rxn_to_fire, dt_macro, da_betabar_dtheta, dNbar_dtheta] = micro_scale(N, input_specs.k, input_specs.fast_rxns, slow_rxns, stoich_fast, input_specs.stoich_react, batch_length_orig, input_specs.num_batches, input_specs.delta);
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
while ind_rec < input_specs.N_record + 1
    
    fprintf(out_file, 'Macroscopic steps executed: %d\n', n_events);               % Eventually write this into an output file
    fprintf(out_file, 'KMC time: %d s\n\n', t);
    
    N_r(ind_rec,:) = N_prev;
    N_int_r(ind_rec,:) = N_int_prev + N_prev * (t_r(ind_rec) - t_prev);
    W_r(ind_rec,:) = W_prev - sum(da_betabar_dtheta) * (t_r(ind_rec) - t_prev);
    dNbar_dtheta_r(:,:,ind_rec) = dNbar_dtheta;
    
    ind_rec = ind_rec + 1;                                          % Increment the recording index
end

elapsed_cpu = cputime-t_cpu_start;
fprintf(out_file, 'Elapsed CPU time: %d seconds\n', elapsed_cpu)

%% Print Data into Output File

fclose('all');
return

fidout = fopen('MSA_output.bin','w');
output_mat = [t_r', N_r, N_int_r, W_r];
fwrite(fidout,output_mat,'double');
fclose(fidout);

fidout = fopen('micro_derivs.bin','w');
fwrite(fidout,dNbar_dtheta_r,'double');
fclose(fidout);

fclose('all');