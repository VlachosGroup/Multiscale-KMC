% STS D-D

function STS_ODE

clc
clear

% Example: Automated TTS D-D simulation

% System specifications
N_0 = [30, 10, 60];                  % initial condition 
t_final = 2;
epsilon = 10^-(3);
k = [1, 1.5, 2, 1, 0.4, 0];
k(1:2) = k(1:2) / epsilon;
S = [1 -1 -1 1 0 0
    0 0 1 -1 -1 1
    -1 1 0 0 1 -1];
[M,R] = size(S);

% Solve the ODE system
[T,Y] = ode15s(@rate_eqns_param,[0, t_final],[N_0,zeros(1,M*R)]);

    function diff = rate_eqns_param(t,z) 
        
        diff = zeros(M + M * R,1);
        N = z(1:M);
        C = z(M+1:end);
        
        [r,dr_dN,dr_dtheta] = rate_eqns(N, k);               % Get rate law information
        diff(1:M) = S * r;                                                           % Get differential changes in species populations
        
        % Indexes the matrix elements of dNdk into a vector for
        % differential equation
        for i = 1:R                   
            diff(i*M+1:(i+1)*M) = S * (dr_dN * C((i-1)*M+1:i*M) + dr_dtheta(:,i));
        end

    end

% Scale the sensitivities by epsilon, this saves us from having to code in
% additional functional dependence on epsilon
Y(:,4:9) = Y(:,4:9) / epsilon;

Y(end,:)'
%% Plot Results
figure
subplot(2,2,1)
plot(T,Y(:,1),T,Y(:,2),T,Y(:,3))
title('Species Populations','FontSize',24)
xlabel('Time (s)','FontSize',24)
ylabel('Population','FontSize',24)
h_legend=legend('A','B','*');
set(h_legend,'FontSize',20);

subplot(2,2,2)
plot(T,Y(:,4),T,Y(:,7),T,Y(:,10),T,Y(:,13),T,Y(:,16),T,Y(:,19))
title('NA Sensitivities','FontSize',24)
xlabel('Time (s)','FontSize',24)
ylabel('Derivative','FontSize',24)
h_legend=legend('k1','k2','k3','k4','k5','k6');
set(h_legend,'FontSize',20);

subplot(2,2,3)
plot(T,Y(:,5),T,Y(:,8),T,Y(:,11),T,Y(:,14),T,Y(:,17),T,Y(:,20))
title('NB Sensitivities','FontSize',24)
xlabel('Time (s)','FontSize',24)
ylabel('Derivative','FontSize',24)
h_legend=legend('k1','k2','k3','k4','k5','k6');
set(h_legend,'FontSize',20);

subplot(2,2,4)
plot(T,Y(:,6),T,Y(:,9),T,Y(:,12),T,Y(:,15),T,Y(:,18),T,Y(:,21))
title('N* Sensitivities','FontSize',24)
xlabel('Time (s)','FontSize',24)
ylabel('Derivative','FontSize',24)
h_legend=legend('k1','k2','k3','k4','k5','k6');
set(h_legend,'FontSize',20);

end