classdef network_input

    properties
        
        % Specify the reaction system - used for all models
        spec_names;                 % names of chemical species
        N_0;                        % initial state
        stoich;                     % stoichiometry matrix
        stoich_react;               % stoichiometric matrix for reactants only
        k;                          % rate constants of elementary reactions
        param_names;                % names of the parameters
        t_final;                    % termination time (s)
        N_record = 1001;            % number of time points to record
        fast_rxns;                  % indices of the fast reactions
        
        % Only for STS ODE and KMC
        eps = 0.1;                % stiffness level, eps << 1
        
        % Only for TTS KMC
        num_batches = 50;        % number of batches to use for microscale steady-state averaging in KMC_TTS
        delta = 0.05;              % accuracy level of microscale averaging, delta << 1
        
    end
    
    methods
        
        function ni = network_input()
            
            ni = ni.AtoB_cat();
            
        end
        
        function ni = plot_species_profiles(ni, t_vec, spec_mat)
            
            [~, n_specs] = size(ni.stoich);
            
            set(0,'defaultlinelinewidth',1.5)
            set(0,'defaultaxeslinewidth',1.5)
            
            % Plot species populations
            figure
            hold on
            for spec = 1:n_specs
                plot(t_vec, spec_mat(:,spec))
            end
            hold off
            box('on')
            ax = gca;
            ax.FontSize = 18;
            xlabel('Time (s)','FontSize',18)
            ylabel('Species count','FontSize',18)
            legend(ni.spec_names);
            legend('boxoff')
            
        end
        
        function ni = AtoB_cat(ni)
            
            ni.spec_names = {'A*','B*','*'};
            ni.N_0 = [0 0 100];
            
            ni.stoich = [1 0 -1
                -1 0 1
                -1 1 0
                1 -1 0
                0 -1 1
                0 1 -1];         % stoichiometric matrix
            ni.k = [4.72, 4.72/6.67, 2, 1, 0.4, 0];                                                         % Rate constants
            
            ni.param_names = {};
            for i = 1:length(ni.k)
                ni.param_names = [ni.param_names ['k_' num2str(i)]];
            end
            
            ni.t_final = 5;
            ni.fast_rxns = [1,2];                                                                       % Stiffness parameter
            
            % Compute reactant matrix
            ni.stoich_react = ni.stoich;
            ni.stoich_react(ni.stoich_react > 0) = 0;
            ni.stoich_react = -ni.stoich_react;
            
        end
        
        function ni = ABC(ni)
            
            ni.spec_names = {'A','B','C'};
            ni.N_0 = [100 0 0];                                                             % Initial state
            
            ni.stoich = [-1 1 0
                1 -1 0
                0 -1 1];         % stoichiometric matrix
            ni.k = [1, 1.5, 2];                                                         % Rate constants%
            ni.param_names = {};
            for i = 1:length(k)
                ni.param_names = [ni.param_names ['k_' num2str(i)]];
            end
            
            ni.t_final = 5;
            ni.fast_rxns = [1,2];
            
            
            % Compute reactant matrix
            ni.stoich_react = ni.stoich_react;
            ni.stoich_react(ni.stoich_react > 0) = 0;
            ni.stoich_react = -ni.stoich_react;
            
        end
        
    end
    
end