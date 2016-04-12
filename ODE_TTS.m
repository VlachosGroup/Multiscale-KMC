% Differential algebraic equation (DAE), two time-scale ODE
% Two time scale (TTS)

function ODE_TTS

clc; clear; fclose('all');

%% System inputs

[spec_names, N_0, stoich, S_react, k, t_final, N_record, fast_rxns, eps, num_batches, delta] = FauxInputRead


N_0 = N_0';     % make it a column vector for convenience
stoich = stoich';
[M,R] = size(stoich);

% Split the stoichiometry matrix into fast and slow components
Ss = stoich;
Ss(:,fast_rxns) = 0;
Sf = stoich - Ss;

% Do some sort of Gaussian elimination procedure on soich_fast to get Tr, Tr_fast,
% and Tr_slow
% Could automate the creation of Ss, Sf, T, Tf, Ts given a specification of
% which reactions are fast/slow

% Transformation matrix for the variables to a new system

Ts = null(Sf','r')';

Tf = [1 0 0]            % How to get this in an automated way? Need to compute Tf so that it is linearly independent from Ts (and is not in the null space of Tf)
T = [Tf; Ts]            % Need to use some sorft of basis set, row echelon form? linear algebra techniques
return

T = [1 0 0; 0 1 0; 1 1 1];                        % Need to be able to compute this transform matrix from the stoich_fast matrix
mf = 1;
ms = M - mf;
Tinv = inv(T);
Tf = T(1:mf,:);
Ts = T(mf+1:end,:);
fse = 2*mf + ms;                        % number of fast-scale equations

% Initial conditions
y_0 = T * N_0;
Cy_0 = zeros(fse*R,1);

%% Solve the ODE system

% Create the mass matrix

num_eqns = M + fse * R;
mass = eye(num_eqns);
mass(1:mf,:) = zeros(mf,num_eqns);                 % Algebraic equations for y_f
% Algebraic equations for C_yf
for j = 1:R
    mass(M+1 + (j-1) * fse : M + (j-1) * fse + 2*mf, :) = zeros(2*mf,num_eqns);
end

% Call ODE solver
options = odeset('Mass',mass);
[T,Y] = ode15s(@rate_eqns_param_TTS,[0 t_final],[y_0;Cy_0],options);

    % Function for ODE solver
    function dz_dt = rate_eqns_param_TTS(t,z)
        
        dz_dt = zeros(num_eqns,1);
        
        % Get rate info
        [rates,dr_dtheta,dr_dN] = rxn_rates(S_react, (Tinv * z(1:M))', k);
        
        % Differentials of y
        dz_dt(1:mf) = Tf * Sf * rates;              % fast
        dz_dt(mf+1:M) = Ts * Ss * rates;            % slow
        
        % Differentials of Cy
        for i = 1:R
            dgdy  = Tf * Sf * dr_dN * Tinv(:,1:mf);
            dz_dt((i-1)*fse+1+M : (i-1)*fse+M+mf) = Tf * Sf * (dr_dN * Tinv * [z((i-1)*fse+1+M :(i-1)*fse+M+mf) + z((i-1)*fse+M+mf+1:(i-1)*fse+M+2*mf); z((i-1)*fse+M+2*mf+1:i*fse+M)]  + dr_dtheta(:,i));  % fast-fast and fast-slow        
            dz_dt((i-1)*fse+M+mf+1 : (i-1)*fse+M+2*mf) = z((i-1)*fse+1+M :(i-1)*fse+M+mf) + inv(dgdy) * Tf * Sf * dr_dtheta(:,i);                                                                                                                   % fast-fast
            dz_dt((i-1)*fse+M+2*mf+1 : i*fse+M) = Ts * Ss * (dr_dN * Tinv * [z((i-1)*fse+1+M :(i-1)*fse+M+mf) + z((i-1)*fse+M+mf+1:(i-1)*fse+M+2*mf); z((i-1)*fse+M+2*mf+1:i*fse+M)] + dr_dtheta(:,i));                                                                                                       % slow
        end
        
        
    end

%% Plot Results - Original Variables

ticklabelfontsize = 20;
axislabelfontsize = 24;



% Transform back to original 
Tinv_aug = blkdiag(Tinv(:,1:mf),Tinv);                                                  % Add additional inverse lines because of the extra variables we have
Tbig = blkdiag(Tinv,Tinv_aug,Tinv_aug,Tinv_aug,Tinv_aug,Tinv_aug,Tinv_aug);             % Need to automate this line so it works for different systems
Yorig = Y * Tbig';

% Return the values of the averages and derivatives at termination time
N_final = Yorig(end,1:3)'
slow_ders = Yorig(end,[7,8,9, 13,14,15, 19,20,21, 25,26,27, 31,32,33]);           
fast_ders = Yorig(end,[4,5,6, 10,11,12, 16,17,18, 22,23,24, 28,29,30]);
ders = slow_ders + fast_ders;
ders = reshape(ders,[3,5])'
% for j1 = 1:5
%     for j2 = 1:3
%         ders(j1,j2) = ders(j1,j2) * k(j1) / Yorig(end,j2);
%     end
% end

% Applies to all graphs
set(0,'defaultlinelinewidth',2)
set(0,'defaultaxeslinewidth',2)
mar = 0.15;

figure
subplot(2,2,1)
plot(T,Yorig(:,1),'k--',T,Yorig(:,2),'k',T,Yorig(:,3),'k-.')
set(gca,'XTickLabel','')
% h_legend=legend('A','B','*');
% h_legend.FontSize = ticklabelfontsize;
set(gca,'Position',[mar 0.5+0.05 0.5-mar 0.5-mar])
ax = gca;
ax.XTick = [1 2 3];
ax.YTick = [ 20 40 60];
ax.FontSize = ticklabelfontsize;
xlim([0 4])
ylim([0 60])
ylabel('Population','FontSize',axislabelfontsize)

%% Complete Derivatives of Species
subplot(2,2,2)
%plot(T,(Yorig(:,4) + Yorig(:,7))*k(1)/Yorig(end,1),T,(Yorig(:,10) + Yorig(:,13))*k(2)/Yorig(end,1),T,(Yorig(:,16) + Yorig(:,19))*k(3)/Yorig(end,1),T,(Yorig(:,22) + Yorig(:,25))*k(4)/Yorig(end,1),T,(Yorig(:,28) + Yorig(:,31))*k(5)/Yorig(end,1),T,(Yorig(:,34) + Yorig(:,37))*k(6)/Yorig(end,1))
%plot(T,(Yorig(:,4) + Yorig(:,7))*k(1)/Yorig(end,1),T,(Yorig(:,10) + Yorig(:,13))*k(2)/Yorig(end,1),T,(Yorig(:,16) + Yorig(:,19))*k(3)/Yorig(end,1),T,(Yorig(:,22) + Yorig(:,25))*k(4)/Yorig(end,1),T,(Yorig(:,28) + Yorig(:,31))*k(5)/Yorig(end,1))
%plot(T,(Yorig(:,4) + Yorig(:,7)),'k',T,(Yorig(:,10) + Yorig(:,13)),'r',T,(Yorig(:,16) + Yorig(:,19)),'m',T,(Yorig(:,22) + Yorig(:,25)),'*g',T,(Yorig(:,28) + Yorig(:,31)),'b')
p2 = plot(T,(Yorig(:,4) + Yorig(:,7)),T,(Yorig(:,10) + Yorig(:,13)),'r',T,(Yorig(:,16) + Yorig(:,19)),'m',T,(Yorig(:,22) + Yorig(:,25)),'g');
p2(1).Color = [0.5 0.5 0.5];
p2(4).Color = [0 1 0];
%title('NA Sensitivities','FontSize',24)
ax = gca;
ax.YAxisLocation = 'right';
ax.FontSize = ticklabelfontsize;
ax.XTick = [1 2 3];
ax.YTick = [0 10 20];
% h_legend=legend('\alpha_1','\alpha_2','\beta_1','\beta_2, \beta_3');
% h_legend.FontSize = ticklabelfontsize;
set(gca,'Position',[0.5 0.5+0.05  0.5-mar 0.5-mar])
xlim([0 4])
ylim([-10 20])
ylabel('Sensitivity','FontSize',axislabelfontsize)

subplot(2,2,3)
%plot(T,(Yorig(:,5) + Yorig(:,8))*k(1)/Yorig(end,2),T,(Yorig(:,11) + Yorig(:,14))*k(2)/Yorig(end,2),T,(Yorig(:,17) + Yorig(:,20))*k(3)/Yorig(end,2),T,(Yorig(:,23) + Yorig(:,26))*k(4)/Yorig(end,2),T,(Yorig(:,29) + Yorig(:,32))*k(5)/Yorig(end,2),T,(Yorig(:,35) + Yorig(:,38))*k(6)/Yorig(end,2))
%plot(T,(Yorig(:,5) + Yorig(:,8))*k(1)/Yorig(end,2),T,(Yorig(:,11) + Yorig(:,14))*k(2)/Yorig(end,2),T,(Yorig(:,17) + Yorig(:,20))*k(3)/Yorig(end,2),T,(Yorig(:,23) + Yorig(:,26))*k(4)/Yorig(end,2),T,(Yorig(:,29) + Yorig(:,32))*k(5)/Yorig(end,2))
%plot(T,(Yorig(:,5) + Yorig(:,8)),'k',T,(Yorig(:,11) + Yorig(:,14)),'r',T,(Yorig(:,17) + Yorig(:,20)),'m',T,(Yorig(:,23) + Yorig(:,26)),'*g',T,(Yorig(:,29) + Yorig(:,32)),'b')
p3 = plot(T,(Yorig(:,5) + Yorig(:,8)),T,(Yorig(:,11) + Yorig(:,14)),'r',T,(Yorig(:,17) + Yorig(:,20)),'m',T,(Yorig(:,23) + Yorig(:,26)),'g');
p3(1).Color = [0.5 0.5 0.5];
set(gca,'Position',[mar mar+0.05 0.5-mar 0.5-mar])
ax = gca;
ax.FontSize = ticklabelfontsize;
ax.XTick = [0 1 2 3 4];
ax.YTick = [-20 -10 0 10 20];
xlim([0 4])
ylim([-20 20])
xlabel('Time (s)')
ylabel('Sensitivity','FontSize',axislabelfontsize)

disp('this is it')
T(end)
Yorig(end,23) + Yorig(end,26)

disp('hi')
subplot(2,2,4)
%plot(T,(Yorig(:,6) + Yorig(:,9))*k(1)/Yorig(end,3),T,(Yorig(:,12) + Yorig(:,15))*k(2)/Yorig(end,3),T,(Yorig(:,18) + Yorig(:,21))*k(3)/Yorig(end,3),T,(Yorig(:,24) + Yorig(:,27))*k(4)/Yorig(end,3),T,(Yorig(:,30) + Yorig(:,33))*k(5)/Yorig(end,3),T,(Yorig(:,36) + Yorig(:,39))*k(6)/Yorig(end,3))
%plot(T,(Yorig(:,6) + Yorig(:,9))*k(1)/Yorig(end,3),T,(Yorig(:,12) + Yorig(:,15))*k(2)/Yorig(end,3),T,(Yorig(:,18) + Yorig(:,21))*k(3)/Yorig(end,3),T,(Yorig(:,24) + Yorig(:,27))*k(4)/Yorig(end,3),T,(Yorig(:,30) + Yorig(:,33))*k(5)/Yorig(end,3))
%plot(T,(Yorig(:,6) + Yorig(:,9)),'k',T,(Yorig(:,12) + Yorig(:,15)),'r',T,(Yorig(:,18) + Yorig(:,21)),'m',T,(Yorig(:,24) + Yorig(:,27)),'*g',T,(Yorig(:,30) + Yorig(:,33)),'b')
p4 = plot(T,(Yorig(:,6) + Yorig(:,9)),'k',T,(Yorig(:,12) + Yorig(:,15)),'r',T,(Yorig(:,18) + Yorig(:,21)),'m',T,(Yorig(:,24) + Yorig(:,27)),'g');
p4(1).Color = [0.5 0.5 0.5];
set(gca,'Position',[0.5 mar+0.05 0.5-mar 0.5-mar])
set(gca,'YAxisLocation','right');
ax = gca;
ax.FontSize = ticklabelfontsize;
ax.XTick = [1 2 3 4];
ax.YTick = [-30 -20 -10 0 10 20];
xlim([0 4])
ylim([-30 20])
xlabel('Time (s)')
ylabel('Sensitivity','FontSize',axislabelfontsize)

set(gcf,'un','n','pos',[0.25,0.25,0.4,0.5]);

%% Fast-scale Derivatives of Species
% figure
% plot(T,Yorig(:,4)*k(1)/Yorig(end,1),T,Yorig(:,10)*k(2)/Yorig(end,1),T,Yorig(:,16)*k(3)/Yorig(end,1),T,Yorig(:,22)*k(4)/Yorig(end,1),T,Yorig(:,28)*k(5)/Yorig(end,1),T,Yorig(:,34)*k(6)/Yorig(end,1))
% title('NA Fast Sensitivities','FontSize',24)
% set(gca,'FontSize',16) %set the font size of everything, including the tick labels
% xlhand = get(gca,'xlabel'); %make a handle for the x axis label
% xlabel('Time (s)') %label the x axis
% set(xlhand,'fontsize',24) 
% ylhand = get(gca,'ylabel'); %make a handle for the y axis label
% ylabel('Derivative') %label the y axis
% set(ylhand,'fontsize',24) %set the font size for the y axis label
% h_legend=legend('k1','k2','k3','k4','k5','k6');
% set(h_legend,'FontSize',20);
% 
% figure
% plot(T,Yorig(:,5)*k(1)/Yorig(end,2),T,Yorig(:,11)*k(2)/Yorig(end,2),T,Yorig(:,17)*k(3)/Yorig(end,2),T,Yorig(:,23)*k(4)/Yorig(end,2),T,Yorig(:,29)*k(5)/Yorig(end,2),T,Yorig(:,35)*k(1)/Yorig(end,2))
% title('NB Fast Sensitivities','FontSize',24)
% set(gca,'FontSize',16) %set the font size of everything, including the tick labels
% xlhand = get(gca,'xlabel'); %make a handle for the x axis label
% xlabel('Time (s)') %label the x axis
% set(xlhand,'fontsize',24) 
% ylhand = get(gca,'ylabel'); %make a handle for the y axis label
% ylabel('Derivative') %label the y axis
% set(ylhand,'fontsize',24) %set the font size for the y axis label
% h_legend=legend('k1','k2','k3','k4','k5','k6');
% set(h_legend,'FontSize',20);
% 
% figure
% plot(T,Yorig(:,6)*k(1)/Yorig(end,3),T,Yorig(:,12)*k(2)/Yorig(end,3),T,Yorig(:,18)*k(3)/Yorig(end,3),T,Yorig(:,24)*k(4)/Yorig(end,3),T,Yorig(:,30)*k(5)/Yorig(end,3),T,Yorig(:,36)*k(6)/Yorig(end,3))
% title('N* Fast Sensitivities','FontSize',24)
% set(gca,'FontSize',16) %set the font size of everything, including the tick labels
% xlhand = get(gca,'xlabel'); %make a handle for the x axis label
% xlabel('Time (s)') %label the x axis
% set(xlhand,'fontsize',24) 
% ylhand = get(gca,'ylabel'); %make a handle for the y axis label
% ylabel('Derivative') %label the y axis
% set(ylhand,'fontsize',24) %set the font size for the y axis label
% h_legend=legend('k1','k2','k3','k4','k5','k6');
% set(h_legend,'FontSize',20);
% 
% % [Yorig(end,7), Yorig(end,13), Yorig(end,19), Yorig(end,25), Yorig(end,31), Yorig(end,37)
% %     Yorig(end,8), Yorig(end,14), Yorig(end,20), Yorig(end,26), Yorig(end,32), Yorig(end,38)
% %     Yorig(end,9), Yorig(end,15), Yorig(end,21), Yorig(end,27), Yorig(end,33), Yorig(end,39)]
% %% Slow-scale Derivatives of Species
% figure
% plot(T,Yorig(:,7)*k(1)/Yorig(end,1),T,Yorig(:,13)*k(2)/Yorig(end,1),T,Yorig(:,19)*k(3)/Yorig(end,1),T,Yorig(:,25)*k(4)/Yorig(end,1),T,Yorig(:,31)*k(5)/Yorig(end,1),T,Yorig(:,37)*k(6)/Yorig(end,1))
% title('NA Slow Sensitivities','FontSize',24)
% set(gca,'FontSize',16) %set the font size of everything, including the tick labels
% xlhand = get(gca,'xlabel'); %make a handle for the x axis label
% xlabel('Time (s)') %label the x axis
% set(xlhand,'fontsize',24) 
% ylhand = get(gca,'ylabel'); %make a handle for the y axis label
% ylabel('Derivative') %label the y axis
% set(ylhand,'fontsize',24) %set the font size for the y axis label
% h_legend=legend('k1','k2','k3','k4','k5','k6');
% set(h_legend,'FontSize',20);
% 
% figure
% plot(T,Yorig(:,8)*k(1)/Yorig(end,2),T,Yorig(:,14)*k(2)/Yorig(end,2),T,Yorig(:,20)*k(3)/Yorig(end,2),T,Yorig(:,26)*k(4)/Yorig(end,2),T,Yorig(:,32)*k(5)/Yorig(end,2),T,Yorig(:,38)*k(6)/Yorig(end,2))
% title('NB Slow Sensitivities','FontSize',24)
% set(gca,'FontSize',16) %set the font size of everything, including the tick labels
% xlhand = get(gca,'xlabel'); %make a handle for the x axis label
% xlabel('Time (s)') %label the x axis
% set(xlhand,'fontsize',24) 
% ylhand = get(gca,'ylabel'); %make a handle for the y axis label
% ylabel('Derivative') %label the y axis
% set(ylhand,'fontsize',24) %set the font size for the y axis label
% h_legend=legend('k1','k2','k3','k4','k5','k6');
% set(h_legend,'FontSize',20);
% 
% figure
% plot(T,Yorig(:,9)*k(1)/Yorig(end,3),T,Yorig(:,15)*k(2)/Yorig(end,3),T,Yorig(:,21)*k(3)/Yorig(end,3),T,Yorig(:,27)*k(4)/Yorig(end,3),T,Yorig(:,33)*k(5)/Yorig(end,3),T,Yorig(:,39)*k(6)/Yorig(end,3))
% title('N* Slow Sensitivities','FontSize',24)
% set(gca,'FontSize',16) %set the font size of everything, including the tick labels
% xlhand = get(gca,'xlabel'); %make a handle for the x axis label
% xlabel('Time (s)') %label the x axis
% set(xlhand,'fontsize',24) 
% ylhand = get(gca,'ylabel'); %make a handle for the y axis label
% ylabel('Derivative') %label the y axis
% set(ylhand,'fontsize',24) %set the font size for the y axis label
% h_legend=legend('k1','k2','k3','k4','k5','k6');
% set(h_legend,'FontSize',20);

end

