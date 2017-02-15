% Reads the output files from MSA
% You also need a version which reads through many files and averages/takes
% sensitivities

% Do for steady-state

clear; clc;

n_params = 5;
n_specs = 3;
N_record = 1001;
cols = 12;                       % t, N1, N2, N3, N1int, N2int, N3int, W1, W2, W3, W4, W5
Y1 = [];
Y2 = [];

% Reads the names of the folders within the directory
basefldr = './example_data';
d = dir(basefldr);
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];         % ignore current and parent directories
havedata = 0;

for file = 1:length(nameFolds)
    
    flname = [basefldr '\' nameFolds{file} '\MSA_output.bin']
    fidread = fopen(flname,'r');
    data1 = fread(fidread,[N_record,cols],'double');
    fclose(fidread);
    
    flname = [basefldr '\' nameFolds{file} '\micro_derivs.bin'];
    fidread = fopen(flname,'r');
    data2 = fread(fidread,[n_params * n_specs * N_record],'double');
    
    fclose('all');
    
    if isempty(data1) || isempty(data2)
        continue
    else
        havedata = havedata + 1;
        data1(1,:) = [];                        % Cut off t = 0 data
        data2 = reshape(data2, [n_params, n_specs, N_record]);
        data2(:,:,1) = [];                      % Cut off t = 0 data
        Y1 = cat(3, Y1, data1);                 % time point, quantitiy, trajectory
        Y2 = cat(4, Y2, data2);                 % , trajectory
    end
    
end

disp([num2str(havedata) ' samples have been processed.'])

%% Compute Trajectory-averaged statistics
Y1(:,5,:) = Y1(:,5,:) ./ Y1(:,1,:);                                             % Divide by time to get time-averages for integral quantities
Y1(:,6,:) = Y1(:,6,:) ./ Y1(:,1,:);
Y1(:,7,:) = Y1(:,7,:) ./ Y1(:,1,:);
%Y1(:,2:7,:) = Y1(:,2:7,:) / N_sites;

traj_avgs = mean(Y1,3);                                                         % Get trajectory averages
T = traj_avgs(:,1);

micro_sens = mean(Y2,4) ;%/ N_sites;                                            % Average and rescale the fast-scale sensitivities
micro_sens = cat(2,micro_sens,micro_sens);
micro_sens = reshape(micro_sens,[30,1000])';


Y1 = permute(Y1,[3,2,1]);

macro_sens = [];
for i = 1:N_record-1
    sens_mat = cov(Y1(:,2:end,i));                                              % At a certain timepoint, take the covariance of everything except time
    macro_sens = cat(3, macro_sens, sens_mat(7:end,1:6));
end

macro_sens = reshape(macro_sens,[30,1000])';
total_sens = micro_sens + macro_sens;

%% Analyze Statistical error in estimates

% Get variances for the macro part of the sensitivities
% t, N1, N2, N3, N1int, N2int, N3int, W1, W2, W3, W4, W5
% 1  2   3   4     5      6      7    8   9   10  11  12

s_1_3s_vars_CLR = zeros(5,3);
s_1_3s_vars_CELR = zeros(5,3);
s_100s_vars_CLR = zeros(5,3);
s_100s_vars_CELR = zeros(5,3);

s_1_3s_CI_CLR = zeros(5,3);
s_1_3s_CI_CELR = zeros(5,3);
s_100s_CI_CLR = zeros(5,3);
s_100s_CI_CELR = zeros(5,3);

N_bs = 1000;

for i = 1:3
    for j = 1:5
        b_test = bootstrp(N_bs,@cov_calc, Y1(:,1+i,13), Y1(:,7+j,13));
        b_test = sort(b_test);
        s_1_3s_vars_CLR(j,i) = var(b_test);
        s_1_3s_CI_CLR(j,i) = (b_test(round(0.975*N_bs)) - b_test(round(0.025*N_bs)))/2; % 95% confidence interval
        
        b_test = bootstrp(N_bs,@cov_calc, Y1(:,4+i,13), Y1(:,7+j,13));
        b_test = sort(b_test);
        s_1_3s_vars_CELR(j,i) = var(b_test);
        s_1_3s_CI_CELR(j,i) = (b_test(round(0.975*N_bs)) - b_test(round(0.025*N_bs)))/2; % 95% confidence interval
        
        b_test = bootstrp(N_bs,@cov_calc, Y1(:,1+i,1000), Y1(:,7+j,1000));
        b_test = sort(b_test);
        s_100s_vars_CLR(j,i) = var(b_test);
        s_100s_CI_CLR(j,i) = (b_test(round(0.975*N_bs)) - b_test(round(0.025*N_bs)))/2; % 95% confidence interval
        
        b_test = bootstrp(N_bs,@cov_calc, Y1(:,4+i,1000), Y1(:,7+j,1000));
        b_test = sort(b_test);
        s_100s_vars_CELR(j,i) = var(b_test);
        s_100s_CI_CELR(j,i) = (b_test(round(0.975*N_bs)) - b_test(round(0.025*N_bs)))/2; % 95% confidence interval
    end
end

% xlswrite('Var.xlsx', s_1_3s_vars_CLR, 's_1_3s_vars_CLR');
% xlswrite('Var.xlsx', s_1_3s_vars_CELR, 's_1_3s_vars_CELR');
% xlswrite('Var.xlsx', s_100s_vars_CLR, 's_100s_vars_CLR');
% xlswrite('Var.xlsx', s_100s_vars_CELR, 's_100s_vars_CELR');

xlswrite('Var.xlsx', s_1_3s_CI_CLR, 's_1_3s_CI_CLR');
xlswrite('Var.xlsx', s_1_3s_CI_CELR, 's_1_3s_CI_CELR');
xlswrite('Var.xlsx', s_100s_CI_CLR, 's_100s_CI_CLR');
xlswrite('Var.xlsx', s_100s_CI_CELR, 's_100s_CI_CELR');

%% Plot

fidout = fopen('T.bin','w');
fwrite(fidout,T,'double');
fclose(fidout);

fidout = fopen('total_sens.bin','w');
fwrite(fidout,total_sens,'double');
fclose(fidout);

set(0,'defaultlinelinewidth',2)
set(0,'defaultaxeslinewidth',2)

ep = 1000;
t_term = T(ep);

% NA Sensitivities - total sensitivity
figure
plot(T(1:ep), total_sens(1:ep,7), 'b', T(1:ep), total_sens(1:ep,22), 'r')
ax = gca;
ax.FontSize = 20;
xlabel('Time (s)','fontsize',24)
ylabel('NA Sensitivities','fontsize',24)  

% NB Sensitivities - total sensitivity
figure
plot(T(1:ep), total_sens(1:ep,8), 'b', T(1:ep), total_sens(1:ep,23), 'r')
ax = gca;
ax.FontSize = 20;
xlabel('Time (s)','fontsize',24)
ylabel('NB Sensitivities','fontsize',24)                                                              % label the y axis

% N* Sensitivities - total sensitivity
figure
plot(T(1:ep), total_sens(1:ep,10), 'b', T(1:ep), total_sens(1:ep,25), 'r')
ax = gca;
ax.FontSize = 20;
h_legend=legend('CLR', 'CELR');
h_legend.FontSize = 20;
xlabel('Time (s)','fontsize',24)
ylabel('N* Sensitivities','fontsize',24)  