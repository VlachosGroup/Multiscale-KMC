% Inputs
% N: state that the macro solver is in, will uniquely identify the fast
% class, acts as the initial state
% k: all rate constants in the problem
% fast reactions: 
% slow reactions: 
% stoich_fast: Only fast reactions will fire, so these are the only ones
% this script needs to worry about

% Outputs
% N average: "state property" for the macro step
% slow propensity averages: 

%(fast and slow, for ALL parameters) will be used for sensitivity analysis
% on the macro scale (two contributions)

function [N_int_avg, N_origin, betabar, slow_rxn_to_fire, dt_macro, da_betabar_dtheta, dNbar_dtheta] = micro_scale(N, k, fast_rxns, slow_rxns, stoich_fast, S_react, batch_length_orig, num_batches, delta)

% This script takes in system information from a macro solver and 

% All reactions are on the same time scale. Use this to analyze the
% effectiveness of the batch stopping criteria.

% The multi-scale sensitivity analysis of (M. Núñez and D.G. Vlachos, J. Chem. Phys. 142 (4), 044108 (2015)) is used.

n_specs = length(N);
n_rxns_fast = length(fast_rxns);
n_rxns_slow = length(slow_rxns);
n_params = n_rxns_fast + n_rxns_slow;           % In general, this does not need to equal the number of reactions

%% Batch parameters
batch_goal = 1;                                             % will converge the first species number, You could also make this any function of the species state, should this also be input?
batch_length = batch_length_orig;                           % We're using our prior knowledge of the relaxation time to know the approximate time-scale
t_final = num_batches * batch_length;                       % increase batch_length until they are uncorrelated, then increase num_batches until things converge
N_record_perbatch_original = 100;                           % sample points, not sure if this is necessary since we will use integral times...
N_record_perbatch = N_record_perbatch_original;
N_record_orig = N_record_perbatch_original * num_batches;
N_record = N_record_perbatch * num_batches;
converged = 0;                                              % At first, it is not converged
corr_tol = 0.05;                                            % How much correlation between adjacent batches should we tolerate?
t_test = 1.734;                                             % DOF = 18, 90% confidence

% Simulation initialization                                                         % An initialization, will be truncated in taking statistics
t = 0;                                                                              % Initial macro time
t_r = linspace(t_final / N_record, t_final, N_record);                              % Recording times for macro sampling, will need to resize this if the simulation adjusts the length of this sampling
ind_rec = 1;                                                                        % Keep track of which time point must be sampled
t_prev = 0;
n_events = 0;

% Species records for macro sampling
N_prev = N;
N_r = zeros(N_record,n_specs);  
N_int = zeros(1,n_specs);
N_int_r = zeros(N_record,n_specs);
N_int_prev = N_int;

% Slow propensities
[a, der] = rxn_rates(S_react, N, k);
a_slow = a(slow_rxns);
a_slow_r = zeros(N_record,n_rxns_slow);
a_slow_int = zeros(1,n_rxns_slow);
a_slow_int_r = zeros(N_record,n_rxns_slow);
a_slow_int_prev = a_slow_int;

% Slow propensity derivatives
da_slow_d_theta = zeros(n_rxns_slow, n_params);
da_slow_d_theta_r = zeros(n_rxns_slow, n_params, N_record);        % will also need an integral value for this
da_slow_d_theta_int = zeros(n_rxns_slow, n_params);  

% Trajectory derivatives
W = zeros(1,n_params);
W_r = zeros(N_record,n_params);
W_prev = W;

%% Stochastic Simulation Loop

% Keep increasing batch size until correlation
% decreases below corr_tol

% To-Do
% Break the data into batches and test correlation
% Also, take the overall average and test the t-test convergence
% Figure out how to adjust t_final, keep the simulation going, and resample
% the data appropriately

growths = 0;

while ~converged
    
    while t < t_final                                          % time horizon
        
        
        % Record the current state as long as time >= t_sample
        while t >= t_r(ind_rec)
            
            % Record the species numbers
            N_r(ind_rec,:) = N_prev;
            N_int_r(ind_rec,:) = N_int_prev + N_prev * (t_r(ind_rec) - t_prev);
            
            % Record the slow propensities
            a_slow_r(ind_rec,:) = a_slow;
            a_slow_int_r(ind_rec,:) = a_slow_int_prev + a_slow * (t_r(ind_rec) - t_prev);
            
            % Record the slow propensity derivatives
            da_slow_d_theta_r(:,:, ind_rec) = da_slow_d_theta;
            
            W_r(ind_rec,:) = W_prev - sum(der(fast_rxns,:)) * (t_r(ind_rec) - t_prev); 
            
            ind_rec = ind_rec + 1;                                                                              % Increment the recording index
        end
        
        % Fire a reaction
        [a, der] = rxn_rates(S_react, N, k);
        a_fast = a(fast_rxns);                                                                                 % For choosing the next reaction
        a_sum = sum(a_fast);
        rxn_ind = min(find(rand(1)<cumsum(a_fast/sum(a_sum))));
        N_prev = N;
        N = N + stoich_fast(rxn_ind,:);
        
        % Time Update
        t_prev = t;
        dt = log(1/rand(1))/a_sum;
        t = t + dt;
        n_events = n_events + 1;
        
        % Record previous values
        N_int_prev = N_int;
        W_prev = W;     
        a_slow_int_prev = a_slow_int;
        
        % Record slow reaction data
        a_slow = a(slow_rxns)';                                                                                 % For averaging for the macro scale, always record it for the current event rather than the previous one
        da_slow_d_theta = der(slow_rxns,:);
        
        % Integral values
        N_int = N_int_prev + dt * N_prev;
        a_slow_int = a_slow_int_prev + dt * a_slow;            
        dW = 1 / a_fast(rxn_ind) * der(fast_rxns(rxn_ind),:) - sum(der(fast_rxns,:)) * dt;        
        W = W + dW;
        
    end
    
    %% Post-process
    
    % Fill in the recording times that were missed
    while ind_rec < N_record + 1        
        
        N_r(ind_rec,:) = N_prev;
        N_int_r(ind_rec,:) = N_int_prev + N_prev * (t_r(ind_rec) - t_prev);       
        
        a_slow_r(ind_rec,:) = a_slow;
        a_slow_int_r(ind_rec,:) = a_slow_int;
        
        W_r(ind_rec,:) = W_prev - sum(der(fast_rxns,:)) * (t_r(ind_rec) - t_prev);            
        
        ind_rec = ind_rec + 1;                                          % Increment the recording index
    end
    
    %% Batch-mean part
    % always throw away the first bin when averaging
    % Indices at the end of each batch
    batch_indices = linspace(N_record_perbatch, N_record, num_batches);
    N_int_perbatch = N_int_r(batch_indices,:);
    batch_avgs = diff(N_int_perbatch) / batch_length;
    f_batches = batch_avgs(:,batch_goal);                        % Put here whatever state function that you are trying to converge, could take this from an external function, could generalize it to converging a vector-valued function   
       
    % Correlation criteria                                      % Was used
    % previously, but typically takes too long to converge in this way
%     mat = cov(f_batches(1:end-1), f_batches(2:end));
%     corr = mat(1,2) / var(f_batches);
    
    % Confidence interval criteria
    err_f = t_test * std(f_batches) / sqrt(num_batches);
    
    % Possible problem: there is a lot of variance in corr
%     if abs(corr) < corr_tol                 % If the batches are uncorrelated, then the simulation is converged
%         converged = 1;
    if err_f / mean(f_batches) < delta                 % If the batches are uncorrelated, then the simulation is converged
        converged = 1;
    else                                    % Otherwise, increase time and keep going. Add on the original time unit
        
        growths = growths + 1;
        
        % Increase final time
        batch_length = batch_length + batch_length_orig;
        t_final = num_batches * batch_length;
        
        % Increase the indexed quantities
        N_record_perbatch = N_record_perbatch + N_record_perbatch_original;
        N_record = N_record_perbatch * num_batches;
        
        % Add space to the recording arrays
        t_r = [t_r linspace(t_final - (N_record_orig - 1) * batch_length_orig / N_record_perbatch_original, t_final, N_record_orig)];
        N_r = [N_r; zeros(N_record_orig, n_specs)];
        N_int_r = [N_int_r; zeros(N_record_orig, n_specs)];
        W_r = [W_r; zeros(N_record_orig, n_params)];  
        a_slow_r = [a_slow_r; zeros(N_record_orig, n_rxns_slow)];
        a_slow_int_r = [a_slow_int_r; zeros(N_record_orig, n_rxns_slow)];
        da_slow_d_theta_r = cat(3, da_slow_d_theta_r, zeros(n_rxns_slow, n_params, N_record_orig));
    end
    
    % if convergence criteria is met, converged = 1
    % if not, increase t_final by some amount and keep going
    
    %mean(f_batches)                    % Use this line to see how the
    %objective function converges towards the analytical value
    
end

%disp(['There were ' num2str(growths) ' growths.'])

%% Steady-state averaging

% Will need to truncate the transient part to make the averaging better.
% Also, use integrals instead of recordings to get better sampling.

% Fast class averaging
N_int_avg = (N_int_r(end,:) - N_int_r(batch_indices(1),:)) / (t_final - t_r(batch_indices(1)));               % Remove batch 1 from the averaging
betabar = (a_slow_int_r(end,:) - a_slow_int_r(batch_indices(1),:)) / (t_final - t_r(batch_indices(1)));               % Remove batch 1 from the averaging                                      
da_betabar_dtheta_C1 = mean(da_slow_d_theta_r(:,:,batch_indices(1):end),3);

% Chose macro reaction and time step
% Choose a State to fire from
a_slow_total_r = sum(a_slow_r,2);                                       % sum of all slow propensities at each step, 
a_slow_total_r(1:N_record_perbatch) = 0;                                % We don't want to choose a state from the first batch because this is unrelaxed
state_to_fire_from = min(find(rand(1)<cumsum(a_slow_total_r/sum(a_slow_total_r))));
N_origin = N_r(state_to_fire_from,:);

% Choose a reaction and time step
a_slow_chosen_state = a_slow_r(state_to_fire_from,:);
slow_rxn_to_fire = min(find(rand(1)<cumsum(a_slow_chosen_state / sum(a_slow_chosen_state))));
dt_macro = log(1/rand(1))/sum(betabar);


%% Sensitivity analysis

% I don't think we can compute the sensitivities with integral values
% unless we sample every event.

% Convert W into Delta_W, use Delta t = batch length
Diff_W = zeros(N_record - N_record_perbatch, n_params);
N_SA = zeros(N_record - N_record_perbatch, n_specs);
a_slow_SA = zeros(N_record - N_record_perbatch, n_rxns_slow);

for i = 1 : N_record - N_record_perbatch
	Diff_W(i,:) = W_r(i + N_record_perbatch,:) - W_r(i,:);
    N_SA(i,:) = N_r(i + N_record_perbatch,:);
    a_slow_SA(i,:) = a_slow_r(i + N_record_perbatch,:);
end

% Create covariance matrices
sens_mat1 = cov([N_SA, Diff_W]);
sens_mat2 = cov([a_slow_SA, Diff_W]);

% Take the parts of the covariance matrices which are the sensitivities we
% want
da_betabar_dtheta_C2 = sens_mat2(1:n_rxns_slow, n_rxns_slow+1:end);                       % Probability rescaling contribution, rows: parameters, cols: slow reactions, should be nonzero for fast parameters and zero for slow parameters
da_betabar_dtheta = da_betabar_dtheta_C1 + da_betabar_dtheta_C2;                        % rows: parameters, cols: slow reactions

dNbar_dtheta = sens_mat1(n_specs+1:end, 1:n_specs);      % Probability rescaling contribution only, rows: parameters, cols: species, this requires covariance = cov(f,Delta_W)

end