# include <iostream>    // using IO functions
# include <string>      // using string
# include <vector>
# include <fstream>
# include <random>
# include <math.h>
# include <sstream>
# include "myHeader.h"
using namespace std;


KMC_traj_TTS :: KMC_traj_TTS() :  KMC_traj(){}            // empty constructor

/*
============= Macroscale KMC simulation =============
*/

void KMC_traj_TTS :: simulate(int rand_seed){       // Execute simulation

srand(rand_seed);      // Set the random seed
KMC_traj::initialize_sim(rand_seed);


N_micro_avg.resize(in_data.n_specs);
for(int k =0; k < in_data.n_specs; k++){
    N_micro_avg[k] = N[k];
}

prop_ders_direct.resize(in_data.n_rxns);
prop_ders_indirect.resize(in_data.n_rxns);
for (int i = 0; i < in_data.n_rxns; i++){
    prop_ders_direct[i].resize(in_data.n_params);
    prop_ders_indirect[i].resize(in_data.n_params);
}


// Vector for microscale sensitivities

micro_scale_sens.resize(in_data.n_specs);
for (int j = 0; j < in_data.n_specs; j++){
    micro_scale_sens[j].resize(in_data.n_params);
}

// Vector for recording microscale sensitivity profiles
micro_scale_sens_profile.resize(in_data.N_record);
for (int i = 0; i < in_data.N_record; i++){
    micro_scale_sens_profile[i].resize(in_data.n_specs);
    for (int j = 0; j < in_data.n_specs; j++){
        micro_scale_sens_profile[i][j].resize(in_data.n_params);
    }
}

// Initializing as 0s here should not be necessary
for(int i = 0; i < in_data.N_record; i++){
    for(int j = 0; j < in_data.n_specs; j++){
        for(int k = 0; k < in_data.n_params; k++){
            micro_scale_sens_profile[i][j][k] = 0;
        }
    }
}


// Resize microscale vectors
dEdth.resize(in_data.n_params);
dEdth_avg.resize(in_data.n_params);
rev_prop_ders.resize(in_data.n_params);
q_cum.resize(in_data.n_micro_steps);
N_candidate.resize(in_data.n_specs);
micro_props.resize(in_data.n_rxns);
N_rec.resize(in_data.n_micro_steps);
for(int i = 0; i < in_data.n_micro_steps; i++){
    N_rec[i].resize(in_data.n_specs);
}
micro_prop_ders.resize(in_data.n_rxns);
for(int i = 0; i < in_data.n_rxns; i++){
    micro_prop_ders[i].resize(in_data.n_params);
}
prop_ders_direct.resize(in_data.n_rxns);
for(int i = 0; i < in_data.n_rxns; i++){
    prop_ders_direct[i].resize(in_data.n_params);
}
prop_ders_indirect.resize(in_data.n_rxns);
for(int i = 0; i < in_data.n_rxns; i++){
    prop_ders_indirect[i].resize(in_data.n_params);
}
slow_props.resize(in_data.n_slow_rxns);
slow_props_cum.resize(in_data.n_slow_rxns);

// Start KMC loop
while(t < in_data.t_final){
    
    /*
    ============================ Compute quantities for current step ============================
    */

    // Use a separete method to handle the microscale averaging
    simulate_micro();
    
    // Add derivatives for each parameter
    for(int k=0; k < in_data.n_params; k++){
        prop_ders_sum[k] = 0;
        for(int i=0; i < in_data.n_rxns; i++){
            prop_ders_sum[k] += prop_ders[i][k];
        }
    }
    
    // If all slow are 0, then exit the while loop
    if(dt < 0 or not std::isfinite(dt)){
        break;
    }
    
    // Record the current state as long as time >= t_sample
    while(t >= in_data.t_rec[ind_rec]){
        record_stats_TTS();
    }
    
    /*
    ============ Fire reaction and update system ==============
    */
    
    // Update species populations
    for(int j = 0; j < in_data.n_specs; j++){
        N[j] += in_data.stoich_mat[rxn_to_fire_ind][j];
    }
    
    // Update clock
    t_prev = t;
    t = t + dt;
    
    // Update trajectory derivatives
    for(int k=0; k < in_data.n_params; k++){
        
        W[k] += prop_ders[rxn_to_fire_ind][k] / props[rxn_to_fire_ind];         // contribution from reaction that fires            
        W[k] -= prop_ders_sum[k] * dt;       // contribution from time step
    }

}

// Fill in recording times that were missed
while(ind_rec < in_data.N_record){
    record_stats_TTS();
}

// Close output files
if(in_data.write_traj_files){
    writer_spec.close();
    writer_SA.close();
}






}

/*
============= Microscale KMC simulation =============
*/

void KMC_traj_TTS :: simulate_micro(){      // Implement with analytical solution for A->B for now

    double r_accept;                      // random mumber between 0 and 1 used to accept or reject the reaction
    double r_timestep;                      // to choose macroscale time step
    
    double E = 0;                    // dimensionless energy of the system
    double del_E;
    double a_fwd;
    double a_rev;
    bool if_accept;

    double Q = 0;       // partition function
    
    int fast_rxn_randi;
    int fast_rxn_to_try;
    int rev_rxn_to_try;
    
    
    // Running counts used for averaging
    for(int j = 0; j < in_data.n_specs; j++){
        N_micro_avg[j] = 0;
    }
    
    for(int j = 0; j < in_data.n_rxns; j++){
        props[j] = 0;
    }
    
    
    for(int i = 0; i < in_data.n_params; i++){
        dEdth[i] = 0;
        dEdth_avg[i] = 0;
    }
    
    
    
    double slow_prop_sum;
    double slow_prop_sum_cum = 0;
    double slow_prop_sum_avg;
    
    /*
    ============= Sensitivity analysis part =============
    */
    
    // Microscale sensitivities
    for(int i = 0; i < in_data.n_specs; i++){
        for(int j = 0; j < in_data.n_params; j++){
            micro_scale_sens[i][j] = 0;
        }
    }
    
    // Derivative of propensities - indirect
    for(int i = 0; i < in_data.n_rxns; i++){
        for(int j = 0; j < in_data.n_params; j++){
            prop_ders_indirect[i][j] = 0;
            prop_ders_direct[i][j] = 0;
        }
    }
    
    
    // Should we scale the number of micro steps by the number of pairs of reversible steps?
    for(int micro_step = 0; micro_step < in_data.n_micro_steps; micro_step++){
        
        // Generate random numbers
        //fast_rxn_randi = micro_step % in_data.n_fast_rxns;      // alternates between which reaction to attempt
        fast_rxn_randi = rand() % in_data.n_fast_rxns;      // picks a random index for a fast reaction
        fast_rxn_to_try = in_data.fast_pairs[fast_rxn_randi][0];
        rev_rxn_to_try = in_data.fast_pairs[fast_rxn_randi][1];
        
        r_accept = ((double) rand() / (RAND_MAX));
        
        
        //cout << "\n" << "species numbers" << endl;
        //for(int i = 0; i < in_data.n_specs; i++){
        //    cout << N[i] << "\t";
        //}
        
        //cout << "\n\n" << "micro_prop_ders" << "\n";
        //for(int i = 0; i < in_data.n_rxns; i++){
        //    for(int j = 0; j < in_data.n_params; j++){
        //        cout << micro_prop_ders << "\t";
        //    }
        //    cout << "\n";
        //}
        
        // Compute reaction propensities
        for(int i = 0; i < in_data.n_rxns; i++){
            
            // Use mass action kinetics equation, a = k * [A]^ma * [B]^mb * ...
            micro_props[i] = in_data.rate_const[i];
            for(int j = 0; j < in_data.n_specs; j++){
                if(in_data.stoich_mat[i][j] < 0){//if this species is a reactant, use it to compute the rate
                    micro_props[i] = micro_props[i] * pow (N[j], -in_data.stoich_mat[i][j]);
                }
            }
            
            // Fill in derivatives for each parameter
            for(int k = 0; k < in_data.n_params; k++){
                if(i==k){
                    
                    micro_prop_ders[i][k] = 1.0;
                    
                    for(int j = 0; j < in_data.n_specs; j++){
                        if(in_data.stoich_mat[i][j] < 0){//if this species is a reactant, use it to compute the rate
                            micro_prop_ders[i][k] = micro_prop_ders[i][k] * pow (N[j], -in_data.stoich_mat[i][j]);
                        }
                    }
                    
                }else{
                    micro_prop_ders[i][k] = 0;
                }
            }
        }
        
        
        
        //cout << "\n\n" << "micro_prop_ders" << "\n";
        //for(int i = 0; i < in_data.n_rxns; i++){
        //    for(int j = 0; j < in_data.n_params; j++){
        //        cout << micro_prop_ders[i][j] << "\t";
        //    }
        //    cout << "\n";
        //}
        
        slow_prop_sum = 0;
        for(int i = 0; i < in_data.n_slow_rxns; i++){
            slow_prop_sum += micro_props[ in_data.slow_rxns[i] ];
        }

        slow_prop_sum_cum += slow_prop_sum;
        
        // Compute forward and reverse propensities
        
        // Update species populations
        a_fwd = micro_props[fast_rxn_to_try];
        for(int j = 0; j < in_data.n_specs; j++){
            N_candidate[j] = N[j] + in_data.stoich_mat[fast_rxn_to_try][j];
        }
        
        a_rev = in_data.rate_const[rev_rxn_to_try];
        for(int j = 0; j < in_data.n_specs; j++){
            if(in_data.stoich_mat[rev_rxn_to_try][j] < 0){//if this species is a reactant, use it to compute the rate
                a_rev = a_rev * pow (N_candidate[j], -in_data.stoich_mat[rev_rxn_to_try][j]);
            }
        }
        
        // Compute derivatives of reverse reaction propensity
        for(int i = 0; i < in_data.n_params; i++){
            rev_prop_ders[i] = 0;
        }
        rev_prop_ders[rev_rxn_to_try] = 1;
        for(int j = 0; j < in_data.n_specs; j++){
            if(in_data.stoich_mat[rev_rxn_to_try][j] < 0){//if this species is a reactant, use it to compute the rate
                rev_prop_ders[rev_rxn_to_try] = rev_prop_ders[rev_rxn_to_try] * pow (N_candidate[j], -in_data.stoich_mat[rev_rxn_to_try][j]);
            }
        }
        
        // Test whether to try to move
        if (a_fwd > 0){
            del_E = -log(a_fwd / a_rev);
            if_accept = (r_accept < a_fwd / a_rev);
        }else{
            del_E = 0;
            if_accept = false;
        }
        
        
        //Record stuff for choosing reaction and state
        for(int spec_ind = 0; spec_ind < in_data.n_specs; spec_ind++){
            N_rec[micro_step][spec_ind] = N[spec_ind];
        }
        
        Q += slow_prop_sum * exp( - E );        // running partition function
        q_cum[micro_step] = Q;                  // sum of little partition functions up until that point
        
        // Record data in cumulative counters
        for(int j = 0; j < in_data.n_specs; j++){
            N_micro_avg[j] += N[j];
        }
        
        for(int j = 0; j < in_data.n_rxns; j++){
            props[j] += micro_props[j];
        }
        
        for(int i = 0; i < in_data.n_rxns; i++){
            for(int j = 0; j < in_data.n_params; j++){
                prop_ders_direct[i][j] += micro_prop_ders[i][j];
            }
        }
        
        // Microscale sensitivities
        for(int i = 0; i < in_data.n_specs; i++){
            for(int j = 0; j < in_data.n_params; j++){
                micro_scale_sens[i][j] += N[i] * (-1 * dEdth[j]);
            }
        }
        
        // Derivative of propensities - indirect
        for(int i = 0; i < in_data.n_rxns; i++){
            for(int j = 0; j < in_data.n_params; j++){
                prop_ders_indirect[i][j] += micro_props[i] * (-1 * dEdth[j]);
            }
        }
        
        for(int i = 0; i < in_data.n_params; i++){
            dEdth_avg[i] += dEdth[i];
        }
        
        
        // Execute the change
        if (if_accept){     // Metropolis criterion
            
            for(int j = 0; j < in_data.n_specs; j++){
                N[j] = N_candidate[j];
            }
            
            E += del_E;         // Update the system energy
            
            // Update energy derivatives
            for(int i = 0; i < in_data.n_params; i++){
                // Need the derivative of the BACKWARDS reaction
                dEdth[i] += -1 / a_fwd * micro_prop_ders[fast_rxn_to_try][i] + 1 / a_rev * rev_prop_ders[i];
            }
            
        } 
        
    }
    
    /*
    ============== Finalize averages ==============
    */
    
    for(int j = 0; j < in_data.n_specs; j++){
        N_micro_avg[j] = N_micro_avg[j] / in_data.n_micro_steps;
    }
    
    for(int j = 0; j < in_data.n_rxns; j++){
        props[j] = props[j] / in_data.n_micro_steps;
    }
    
    // Set propensities of fast reactions to zero
    for(int j = 0; j < in_data.n_fast_rxns; j++){
        props[ in_data.fast_pairs[j][0] ] = 0;
    }
    
    // Direct microscale propensity derivatives
    for(int i = 0; i < in_data.n_rxns; i++){
        for(int j = 0; j < in_data.n_params; j++){
            prop_ders_direct[i][j] = prop_ders_direct[i][j] / in_data.n_micro_steps;
        }
    }
    
    // Set slow propensity values equal to zero for fast parameters
    for(int i = 0; i < in_data.n_fast_rxns; i++){
        for(int j = 0; j < in_data.n_params; j++){
            prop_ders_direct[ in_data.fast_pairs[i][0] ][j] = 0;
            prop_ders_indirect[ in_data.fast_pairs[i][0] ][j] = 0;
        }
    }
    
    slow_prop_sum_avg = slow_prop_sum_cum / in_data.n_micro_steps;
    
    for(int i = 0; i < in_data.n_params; i++){
        dEdth_avg[i] = dEdth_avg[i] / in_data.n_micro_steps;
    }
    
    
    // Microscale sensitivities
    for(int i = 0; i < in_data.n_specs; i++){
        for(int j = 0; j < in_data.n_params; j++){
            micro_scale_sens[i][j] = micro_scale_sens[i][j] / in_data.n_micro_steps - N_micro_avg[i] * -1 * dEdth_avg[j];
        }
    }
    
    // Derivative of propensities - indirect
    for(int i = 0; i < in_data.n_rxns; i++){
        for(int j = 0; j < in_data.n_params; j++){
            prop_ders_indirect[i][j] = prop_ders_indirect[i][j] / in_data.n_micro_steps - props[i] * -1 * dEdth_avg[j];
        }
    }
    
    // Derivative of propensities - total (sum of direct and indirect)
    for (int i = 0; i < in_data.n_rxns; i++){
        for (int j = 0; j < in_data.n_params; j++){
            prop_ders[i][j] = prop_ders_direct[i][j] + prop_ders_indirect[i][j];
        }
    }
       
    
    /*
    ============== Choose things for slow reaction ==============
    */

    // Choose macro time step
    r_timestep = ((double) rand() / (RAND_MAX));
    
    if(slow_prop_sum_avg > 0){
        dt = - log(r_timestep) / slow_prop_sum_avg;
    }else{
        dt = -1;
    }
    
    
    // Choose state to fire from
    
    for(int micro_step = 0; micro_step < in_data.n_micro_steps; micro_step++){
        q_cum[micro_step] = q_cum[micro_step] / Q;      // convert to a fraction
    }
    
    double r_state_choose = ((double) rand() / (RAND_MAX));
    int state_to_choose = 0;
    while (  q_cum[state_to_choose] < r_state_choose){
        state_to_choose++;
    }
    
    for(int i = 0; i < in_data.n_specs; i++){
        N[i] = N_rec[state_to_choose][i];
    }
    
    
    // Choose slow reaction to fire
    
    // Compute slow propensities for the state to fire from
    // Use mass action kinetics equation, a = k * [A]^ma * [B]^mb * ...

    for(int k = 0; k < in_data.n_slow_rxns; k++){
        
        slow_props_cum[k] = 0;
        
        int i2 = in_data.slow_rxns[k];
        slow_props[k] = in_data.rate_const[i2];
        for(int j = 0; j < in_data.n_specs; j++){
            if(in_data.stoich_mat[i2][j] < 0){       //if this species is a reactant, use it to compute the rate
                slow_props[k] = micro_props[k] * pow (N[j], -in_data.stoich_mat[i2][j]);
            }
        }
        
        for(int kleq = 0; kleq <= k; kleq++){
            slow_props_cum[k] += slow_props[kleq];
        }
    }
    
    for(int i = 0; i < in_data.n_slow_rxns; i++){
        slow_props_cum[i] = slow_props_cum[i] / slow_props_cum[in_data.n_slow_rxns - 1];
    }

    double r_slow_rxn = ((double) rand() / (RAND_MAX));
    int slow_ind = 0;
    while ( slow_props_cum[slow_ind] < r_slow_rxn ){
        slow_ind++;
    }
    
    rxn_to_fire_ind = in_data.slow_rxns[slow_ind];
      
    
    /*
    ================== Debug by printing out calculated values ==================
    */
    
    bool print_micro_debug = false;
    
    if (print_micro_debug){
        
        cout << endl;
        cout << " Time (s) " << endl;
        cout << t << endl;
        cout << endl;
        
        
        cout << " Species averages " << endl;
        for(int i = 0; i < in_data.n_specs; i++){
            cout << N_micro_avg[i] << "\t";
        }
        cout << endl;
        
        cout << endl;
        cout << " Final state " << endl;
        for(int i = 0; i < in_data.n_specs; i++){
            cout << N[i] << "\t";
        }
        cout << endl;
        
        cout << endl;
        cout << " Reaction to fire " << endl;
        cout << rxn_to_fire_ind << endl;
        
        cout << endl;
        cout << " Average propensities " << endl;
        for(int i = 0; i < in_data.n_rxns; i++){
            cout << props[i] << "\t";
        }
        cout << endl;
        
        cout << endl;
        cout << " Direct propensity averages " << endl;
        for(int i = 0; i < in_data.n_rxns; i++){
            for(int j = 0; j < in_data.n_params; j++){
                cout << prop_ders_direct[i][j] << "\t";
            }
            cout << endl;
        }
        
        cout << endl;
        cout << " Indirect propensity averages " << endl;
        for(int i = 0; i < in_data.n_rxns; i++){
            for(int j = 0; j < in_data.n_params; j++){
                cout << prop_ders_indirect[i][j] << "\t";
            }
            cout << endl;
        }
    
    }
    
}


/*
========= Record the current state =========
*/
void KMC_traj_TTS :: record_stats_TTS(){
    
    t_trunc = in_data.t_rec[ind_rec] - t_prev;
        
    // Record data in output files
    if(in_data.write_traj_files){
        
        // Record the species numbers
        writer_spec << in_data.t_rec[ind_rec] << "\t";
        for(int j = 0; j < in_data.n_specs; j++){
            writer_spec << N_micro_avg[j] << "\t";
        }
        writer_spec << endl;
    
        // Record trajectory derivatives
        writer_SA << in_data.t_rec[ind_rec] << "\t";
        for(int k = 0; k < in_data.n_params; k++){
            writer_SA << W[k] - prop_ders_sum[k] * t_trunc << "\t";
        }
        writer_SA << endl;
    }
    
    // Record data in variables
    
    // Record species populations
    for(int j = 0; j < in_data.n_specs; j++){
        spec_profile[ind_rec][j] = N_micro_avg[j];
        }
    
    // Record trajectory derivatives
    for(int k = 0; k < in_data.n_params; k++){
        traj_deriv_profile[ind_rec][k] = W[k] - prop_ders_sum[k] * t_trunc;
        }
        
    // Record microscale derivatives
    for(int j = 0; j < in_data.n_specs; j++){
        for(int k = 0; k < in_data.n_params; k++){
            micro_scale_sens_profile[ind_rec][j][k] = micro_scale_sens[j][k];
        }
    }   
        
    
    ind_rec += 1;
    
}

double KMC_traj_TTS :: get_micro_scale_sens_profile(int i, int j, int k){
    return micro_scale_sens_profile[i][j][k];
}