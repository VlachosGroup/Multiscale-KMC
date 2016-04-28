%% Plots

% Read in analysis data

set(0,'defaultlinelinewidth',1.5)
set(0,'defaultaxeslinewidth',2)

% NA Sensitivities - total sensitivity
figure
plot(T, total_sens(:,7), 'b', T, total_sens(:,22), 'r')
ax = gca;
ax.FontSize = 20;
xlabel('time (s)','fontsize',24)
ylabel('NA Sensitivities','fontsize',24)  

% NB Sensitivities - total sensitivity
figure
plot(T, total_sens(:,8), 'b', T, total_sens(:,23), 'r')
ax = gca;
ax.FontSize = 20;
xlabel('time (s)','fontsize',24)
ylabel('NB Sensitivities','fontsize',24)

% N* Sensitivities - total sensitivity
figure
plot(T, total_sens(:,10), 'b', T, total_sens(:,25), 'r')
ax = gca;
ax.FontSize = 20;
h_legend=legend('CLR', 'CELR');
h_legend.FontSize = 20;
xlabel('time (s)','fontsize',24)
ylabel('N* Sensitivities','fontsize',24)  