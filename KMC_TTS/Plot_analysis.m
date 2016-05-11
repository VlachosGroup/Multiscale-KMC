% Read in analysis data
clear; clc; fclose('all');
basefldr = textread('Directory.txt','%c')';
addpath('../Network');
[spec_names, N_0, stoich, S_react, k, param_names, t_final, N_record, fast_rxns, eps, num_batches, delta] = FauxInputRead2;
[n_rxns, n_specs] = size(stoich);

fidout = fopen([basefldr '/T.bin'],'r');
T = fread(fidout,[N_record,1],'double');
fclose(fidout);

fidout = fopen([basefldr '/spec_avgs.bin'],'r');
spec_avgs = fread(fidout,[N_record,2*n_specs],'double');
fclose(fidout);

fidout = fopen([basefldr '/total_sens.bin'],'r');
total_sens = fread(fidout,[N_record, 2 * n_specs * n_rxns],'double');
fclose(fidout);

set(0,'defaultlinelinewidth',1.5)
set(0,'defaultaxeslinewidth',2)

% Plot species populations
figure
hold on
for spec = 1:n_specs
    plot(T,spec_avgs(:,spec))
end
hold off
box('on')
ax = gca;
ax.FontSize = 18;
xlabel('time (s)','FontSize',18)
ylabel('spec. pop.','FontSize',18)
legend(spec_names);

% Plot sensitivities for each species
for spec = 1:n_specs
    
    figure
    hold on
    for param = 1:n_rxns
        plot(T, total_sens(:, param + n_rxns * (spec-1) ))
    end
    
    hold off
    box('on')
    xlabel('time (s)','FontSize',18)
    ylabel([spec_names{spec} ' sensitivities'],'FontSize',18)
    ax = gca;
    ax.FontSize = 18;
    legend(param_names);
end