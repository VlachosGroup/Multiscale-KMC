% Post-processing for a single run
% Reads the output files from MSA

clear; clc;

n_params = 5;
n_specs = 3;
N_record = 1001;
cols = 12;                       % t, N1, N2, N3, N1int, N2int, N3int, W1, W2, W3, W4, W5

fidread = fopen('MSA_output.bin','r');
Y = fread(fidread,[N_record,cols],'double');
fclose(fidread);
Y(1,:) = [];        % Cut off t = 0 data

fidread = fopen('micro_derivs.bin','r');
Y2 = fread(fidread,[n_params * n_specs * N_record],'double');
Y2 = reshape(Y2, [n_params, n_specs, N_record]);
fclose('all');
Y2(:,:,1) = [];        % Cut off t = 0 data

%% Draw Graphs

% Species Populations
figure
plot(Y(:,1), Y(:,2), Y(:,1), Y(:,3), Y(:,1), Y(:,4)) 
title('Species Populations','FontSize',24)
set(gca,'FontSize',16)                                                          % set the font size of everything, including the tick labels
xlhand = get(gca,'xlabel');                                                     % make a handle for the x axis label
xlabel('Time (s)')                                                              % label the x axis
set(xlhand,'fontsize',24) 
ylhand = get(gca,'ylabel');                                                     % make a handle for the y axis label
ylabel('Species Populations')                                                   % label the y axis
set(ylhand,'fontsize',24)                                                       % set the font size for the y axis label
h_legend=legend('NA', 'NB', '*');
set(h_legend,'FontSize',20);

% Species Populations Ergodic Averages
figure
plot(Y(:,1), Y(:,5)./Y(:,1), Y(:,1), Y(:,6)./Y(:,1), Y(:,1), Y(:,7)./Y(:,1)) 
title('Species Populations Integrals','FontSize',24)
set(gca,'FontSize',16)                                                          % set the font size of everything, including the tick labels
xlhand = get(gca,'xlabel');                                                     % make a handle for the x axis label
xlabel('Time (s)')                                                              % label the x axis
set(xlhand,'fontsize',24) 
ylhand = get(gca,'ylabel');                                                     % make a handle for the y axis label
ylabel('Species Populations')                                                   % label the y axis
set(ylhand,'fontsize',24)                                                       % set the font size for the y axis label
h_legend=legend('NA', 'NB', '*');
set(h_legend,'FontSize',20);

% Trajectory Derivatives
figure
plot(Y(:,1), Y(:,8), Y(:,1), Y(:,9), Y(:,1), Y(:,10), Y(:,1), Y(:,11), Y(:,1), Y(:,12)) 
title('Trajectory Derivatives','FontSize',24)
set(gca,'FontSize',16)                                                          % set the font size of everything, including the tick labels
xlhand = get(gca,'xlabel');                                                     % make a handle for the x axis label
xlabel('Time (s)')                                                              % label the x axis
set(xlhand,'fontsize',24) 
ylhand = get(gca,'ylabel');                                                     % make a handle for the y axis label
ylabel('Trajectory Derivatives (W)')                                                   % label the y axis
set(ylhand,'fontsize',24)                                                       % set the font size for the y axis label
h_legend=legend('k1', 'k2', 'k3', 'k4', 'k5');
set(h_legend,'FontSize',20);

% Micro Derivatives
figure
plot(Y(:,1), reshape(Y2(1,1,:),[N_record-1,1]), Y(:,1), reshape(Y2(1,1,:),[N_record-1,1]))
title('Micro Derivatives','FontSize',24)
set(gca,'FontSize',16)                                                          % set the font size of everything, including the tick labels
xlhand = get(gca,'xlabel');                                                     % make a handle for the x axis label
xlabel('Time (s)')                                                              % label the x axis
set(xlhand,'fontsize',24) 
ylhand = get(gca,'ylabel');                                                     % make a handle for the y axis label
ylabel('Species Populations')                                                   % label the y axis
set(ylhand,'fontsize',24)                                                       % set the font size for the y axis label
h_legend=legend('dNAdk1 ', 'dN*dk2');
set(h_legend,'FontSize',20);