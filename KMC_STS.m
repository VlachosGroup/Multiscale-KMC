% single time-scale stochastic simulation

clc; clear;

% System info
input_specs = network_input();

input_specs.k(input_specs.fast_rxns) = input_specs.k(input_specs.fast_rxns) / input_specs.eps;
[n_params, n_specs] = size(input_specs.stoich);
N_r = zeros(input_specs.N_record,n_specs); 
N_int = zeros(1,n_specs);
N_int_r = zeros(input_specs.N_record,n_specs);
N_int_prev = N_int;
t_r = linspace(0, input_specs.t_final, input_specs.N_record);                                               % Recording times, remains constant

%% Stochastic Simulation Loop

% Simulation initialization
N = input_specs.N_0;
N_prev = N;                                                                         % An initialization, will be truncated in taking statistics
t = 0;                                                                              % Initial macro time
N_int = zeros(1,n_specs);
ind_rec = 1;                                                                        % Keep track of which time point must be sampled
t_prev = 0;
n_events = 0;

% Sensitivity analysis Parameters
da_dtheta = zeros(n_params,n_params);
W = zeros(1,n_params);
W_r = zeros(input_specs.N_record,n_params);
W_prev = W;

while t < input_specs.t_final                                                                   % (macro) Termination time controls the sampling
    
    % Record the current state as long as time >= t_sample
    while t >= t_r(ind_rec)                                                     % If you record after you compute the next step, but before you update data, then you won't need all the prev variables
        
        disp(['Event # ' num2str(n_events)])
        disp(['Time: ' num2str(t) ' s'])
        
        % Record the species numbers
        N_r(ind_rec,:) = N_prev;
        N_int_r(ind_rec,:) = N_int_prev + N_prev * (t_r(ind_rec) - t_prev);
        W_r(ind_rec,:) = W_prev - sum(da_dtheta) * (t_r(ind_rec) - t_prev);
        
        ind_rec = ind_rec + 1;                                                                              % Increment the recording index
    end
        
    [a, da_dtheta, ~] = rxn_rates(input_specs.stoich_react, N, input_specs.k);
    a_sum = sum(a);
    rxn_to_fire = min(find(rand(1)<cumsum(a/sum(a_sum))));     
    N_prev = N;
    N = N + input_specs.stoich(rxn_to_fire,:);                                 % which reaction will fire? And from which state?
    
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
while ind_rec < input_specs.N_record + 1
    
    N_r(ind_rec,:) = N_prev;
    N_int_r(ind_rec,:) = N_int_prev + N_prev * (t_r(ind_rec) - t_prev);
    W_r(ind_rec,:) = W_prev - sum(da_dtheta) * (t_r(ind_rec) - t_prev);
    
    ind_rec = ind_rec + 1;                                          % Increment the recording index
end

%% Plot species profile

% Plot species profiles
input_specs.plot_species_profiles(t_r',N_r);

%% Print Data into Output File

fidout = fopen('KMC_STS_out.bin','w');
output_mat = [t_r', N_r, N_int_r, W_r];
fwrite(fidout,output_mat,'double');
fclose(fidout);