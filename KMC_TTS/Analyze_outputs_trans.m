% Reads the output files from MSA
% You also need a version which reads through many files and averages/takes
% sensitivities

function Analyze_outputs_trans

clear; clc; fclose('all');
addpath('../Network');
basefldr = textread('Directory.txt','%c')';
[spec_names, N_0, stoich, S_react, k, param_names, t_final, N_record, fast_rxns, eps, num_batches, delta] = FauxInputRead2;

[n_params, n_specs] = size(stoich);
Y1 = [];
Y2 = [];

% Reads the names of the folders within the directory
d = dir(basefldr);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];         % ignore current and parent directories
havedata = 0;

for file = 1:length(nameFolds)
    
    [data1, data2] = Read_KMCTTS_output([basefldr '/' nameFolds{file}], N_record, n_params, n_specs);
    
    fclose('all');
    
    if isempty(data1) || isempty(data2)
        continue
    else
        havedata = havedata + 1;
        data2 = reshape(data2, [n_params, n_specs, N_record]);
        Y1 = cat(3, Y1, data1);                 % time point, quantitiy, trajectory
        Y2 = cat(4, Y2, data2);                 % , trajectory
    end
    
end

fprintf([num2str(havedata) ' samples have been processed.\n'])
if havedata < 3
    error('Too few data samples to do statistics')
end

%% Compute Trajectory-averaged statistics

% Compute ergodic averages
Y1(1, 2 + n_specs : 1 + 2 *n_specs, :) = Y1(1, 2 : 1 + n_specs, :);
for i = 1:n_specs
    Y1(2:end,1 + n_specs + i,:) = Y1(2:end,1 + n_specs + i,:) ./ Y1(2:end,1,:);
end

% Compute averages across trajectories
traj_avgs = mean(Y1,3);                                                         % Get trajectory averages
T = traj_avgs(:,1);     % time points

% Compute microscopic derivatives by averaging across trajectories
micro_sens = mean(Y2,4) ;
micro_sens = cat(2,micro_sens,micro_sens);          % I do not remember why this line is here...
micro_sens = reshape(micro_sens,[2 * n_specs * n_params, N_record])';

% Compute macroscopic derivatives by taking the covariance of species
% numbers and macro-trajectory derivatives
Y1 = permute(Y1,[3,2,1]);
macro_sens = [];
for i = 1:N_record
    sens_mat = cov(Y1(:,2:end,i));                                              % At a certain timepoint, take the covariance of everything except time
    macro_sens = cat(3, macro_sens, sens_mat(2*n_specs+1 : end, 1 : 2*n_specs));
end
macro_sens = reshape(macro_sens,[2 * n_specs * n_params, N_record])';

% Total sensitivity is a sum of micro and macro contributions
total_sens = micro_sens + macro_sens;

%% Analyze Statistical error in estimates

% Get variances for the macro part of the sensitivities
% t, N1, N2, N3, N1int, N2int, N3int, W1, W2, W3, W4, W5
% 1  2   3   4     5      6      7    8   9   10  11  12

% s_1_3s_vars_CLR = zeros(5,3);
% s_1_3s_vars_CELR = zeros(5,3);
% s_100s_vars_CLR = zeros(5,3);
% s_100s_vars_CELR = zeros(5,3);
%
% s_1_3s_CI_CLR = zeros(5,3);
% s_1_3s_CI_CELR = zeros(5,3);
% s_100s_CI_CLR = zeros(5,3);
% s_100s_CI_CELR = zeros(5,3);
%
% N_bs = 1000;
%
% for i = 1:n_specs
%     for j = 1:n_params
%         b_test = bootstrp(N_bs,@cov_calc, Y1(:,1+i,13), Y1(:,7+j,13));
%         b_test = sort(b_test);
%         s_1_3s_vars_CLR(j,i) = var(b_test);
%         s_1_3s_CI_CLR(j,i) = (b_test(round(0.975*N_bs)) - b_test(round(0.025*N_bs)))/2; % 95% confidence interval
%
%         b_test = bootstrp(N_bs,@cov_calc, Y1(:,4+i,13), Y1(:,7+j,13));
%         b_test = sort(b_test);
%         s_1_3s_vars_CELR(j,i) = var(b_test);
%         s_1_3s_CI_CELR(j,i) = (b_test(round(0.975*N_bs)) - b_test(round(0.025*N_bs)))/2; % 95% confidence interval
%
%         b_test = bootstrp(N_bs,@cov_calc, Y1(:,1+i,1000), Y1(:,7+j,1000));
%         b_test = sort(b_test);
%         s_100s_vars_CLR(j,i) = var(b_test);
%         s_100s_CI_CLR(j,i) = (b_test(round(0.975*N_bs)) - b_test(round(0.025*N_bs)))/2; % 95% confidence interval
%
%         b_test = bootstrp(N_bs,@cov_calc, Y1(:,4+i,1000), Y1(:,7+j,1000));
%         b_test = sort(b_test);
%         s_100s_vars_CELR(j,i) = var(b_test);
%         s_100s_CI_CELR(j,i) = (b_test(round(0.975*N_bs)) - b_test(round(0.025*N_bs)))/2; % 95% confidence interval
%     end
% end

%% Write output files from the analysis
fidout = fopen([basefldr '/T.bin'],'w');
fwrite(fidout,T,'double');
fclose(fidout);

fidout = fopen([basefldr '/spec_avgs.bin'],'w');
fwrite(fidout,traj_avgs(:, 2 : 1+2*n_specs),'double');
fclose(fidout);

fidout = fopen([basefldr '/total_sens.bin'],'w');
fwrite(fidout,total_sens,'double');
fclose(fidout);

end

% Reads the output files from MSA

function [Y, Y2] = Read_KMCTTS_output(dir, N_record, n_params, n_specs)

cols = 1 + 2 * n_specs + n_params;                       % t, species, integral species, trajectory derivatives

fidread = fopen([dir '/MSA_output.bin'],'r');
Y = fread(fidread,[N_record,cols],'double');
fclose(fidread);

fidread = fopen([dir '/micro_derivs.bin'],'r');
Y2 = fread(fidread, n_params * n_specs * N_record, 'double');       % parameters * species * time points
Y2 = reshape(Y2, [n_params, n_specs, N_record]);
fclose('all');

end