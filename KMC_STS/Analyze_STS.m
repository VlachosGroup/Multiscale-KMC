% Reads the output files from MSA
% You also need a version which reads through many files and averages/takes
% sensitivities

clear; clc;

n_params = 5;
n_specs = 3;
N_record = 1001;
cols = 12;                       % t, N1, N2, N3, N1int, N2int, N3int, W1, W2, W3, W4, W5

fidread = fopen('KMC_STS_output.bin','r');
Y = fread(fidread,[N_record,cols],'double');
fclose(fidread);
Y(1,:) = [];        % Cut off t = 0 data


%% Draw Graphs

% Species Populations
figure
plot(Y(:,1), Y(:,2), Y(:,1), Y(:,3), Y(:,1), Y(:,4)) 
set(gca,'FontSize',16)                                                          % set the font size of everything, including the tick labels
xlhand = get(gca,'xlabel');                                                     % make a handle for the x axis label
%xlabel('Time (s)')                                                              % label the x axis
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
ylabel('species pop.')                                                   % label the y axis
set(ylhand,'fontsize',24)                                                       % set the font size for the y axis label
h_legend=legend('NA', 'NB', '*');
set(h_legend,'FontSize',20);

% Trajectory Derivatives
figure
%plot(Y(:,1), Y(:,8), Y(:,1), Y(:,9), Y(:,1), Y(:,10), Y(:,1), Y(:,11), Y(:,1), Y(:,12))
plot(Y(:,1), Y(:,8)*10, Y(:,1), Y(:,9)*10, Y(:,1), Y(:,10), Y(:,1), Y(:,11), Y(:,1), Y(:,12))
set(gca,'FontSize',16)                                                          % set the font size of everything, including the tick labels
xlhand = get(gca,'xlabel');                                                     % make a handle for the x axis label
xlabel('time (s)')                                                              % label the x axis
set(xlhand,'fontsize',24) 
ylhand = get(gca,'ylabel');                                                     % make a handle for the y axis label
ylabel('W_k')                                                   % label the y axis
set(ylhand,'fontsize',24)                                                       % set the font size for the y axis label
h_legend=legend('k1', 'k2', 'k3', 'k4', 'k5');
set(h_legend,'FontSize',20);