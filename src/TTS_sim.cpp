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

void KMC_traj_TTS :: simulate_TTS(int rand_seed){       // Execute simulation

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
    
// Some additional non-class varaiables...
double asum;

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
    
    // Sum the slow propensities
    asum = 0;
    for(int i = 0; i < in_data.n_rxns; i++){
        asum += props[i];
    }
    
    // If all propensities are 0, then exit the while loop
    if(asum == 0){
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
    
    int fast_rxn_randi;
    int fast_rxn_to_try;
    int rev_rxn_to_try;
    int N_candidate[in_data.n_specs];
    
    double micro_props[in_data.n_rxns];
    double micro_prop_ders[in_data.n_rxns][in_data.n_params];
    
    // Running counts used for averaging
    double N_cum[in_data.n_specs];
    for(int j = 0; j < in_data.n_specs; j++){
        N_cum[j] = 0;
    }
    
    double props_cum[in_data.n_rxns];
    for(int j = 0; j < in_data.n_rxns; j++){
        props_cum[j] = 0;
    }
    
    double micro_prop_ders_cum[in_data.n_rxns][in_data.n_params];
    for(int i = 0; i < in_data.n_rxns; i++){
        for(int j = 0; j < in_data.n_params; j++){
            micro_prop_ders_cum[i][j] = 0;
        }
    }
    
    double prop_ders_direct[in_data.n_rxns][in_data.n_params];
    double prop_ders_indirect[in_data.n_rxns][in_data.n_params];
    
    // Should we scale the number of micro steps by the number of pairs of reversible steps?
    for(int micro_step = 0; micro_step < in_data.n_micro_steps; micro_step++){
        
        // Generate random numbers
        fast_rxn_randi = rand() % in_data.n_fast_rxns;      // picks a random index for a fast reaction
        fast_rxn_to_try = in_data.fast_pairs[fast_rxn_randi][0];
        rev_rxn_to_try = in_data.fast_pairs[fast_rxn_randi][1];
        
        r_accept = ((double) rand() / (RAND_MAX));
        
        
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
                            micro_prop_ders[i][k] = prop_ders[i][k] * pow (N[j], -in_data.stoich_mat[i][j]);
                        }
                    }
                    
                }else{
                    micro_prop_ders[i][k] = 0;
                }
            }
        }
        
        // Compute forward and reverse propensities
        a_fwd = micro_props[fast_rxn_to_try];
        for(int j = 0; j < in_data.n_specs; j++){
            N_candidate[j] = N[j];
        }
        
        // Update species populations
        for(int j = 0; j < in_data.n_specs; j++){
            N_candidate[j] += in_data.stoich_mat[fast_rxn_to_try][j];
        }
        
        a_rev = in_data.rate_const[fast_rxn_to_try];
        for(int j = 0; j < in_data.n_specs; j++){
            if(in_data.stoich_mat[fast_rxn_to_try][j] < 0){//if this species is a reactant, use it to compute the rate
                a_rev = a_rev * pow (N_candidate[j], -in_data.stoich_mat[fast_rxn_to_try][j]);
            }
        }
        
        // Test whether to try to move
        del_E = -log(a_fwd / a_rev);
        if_accept = (r_accept < a_fwd / a_rev);
        
        // Execute the change
        if (if_accept){     // Metropolis criterion
            
            for(int j = 0; j < in_data.n_specs; j++){
                N[j] = N_candidate[j];
            }
            
            E += del_E;         // Update the system energy
            
        }  
        
        // Record data in cumulative counters
        for(int j = 0; j < in_data.n_specs; j++){
            N_cum[j] += N[j];
        }
        
        for(int j = 0; j < in_data.n_rxns; j++){
            props_cum[j] += micro_props[j];
        }
        
        for(int i = 0; i < in_data.n_rxns; i++){
            for(int j = 0; j < in_data.n_params; j++){
                micro_prop_ders_cum[i][j] += micro_prop_ders[i][j];
            }
        }
        
    }
    
    // Finalize averages
    for(int j = 0; j < in_data.n_specs; j++){
        N_micro_avg[j] = N_cum[j] / in_data.n_micro_steps;
    }
    
    for(int j = 0; j < in_data.n_rxns; j++){
        props[j] = props_cum[j] / in_data.n_micro_steps;
    }
    
    // Set propensities of fast reactions to zero
    for(int j = 0; j < in_data.n_fast_rxns; j++){
        props[ in_data.fast_pairs[j][0] ] = 0;
    }
    
    // Direct microscale propensity derivatives
    for(int i = 0; i < in_data.n_rxns; i++){
        for(int j = 0; j < in_data.n_params; j++){
            prop_ders_direct[i][j] += micro_prop_ders_cum[i][j] / in_data.n_micro_steps;
        }
    }
    
    for(int i = 0; i < in_data.n_fast_rxns; i++){
        for(int j = 0; j < in_data.n_params; j++){
            prop_ders_direct[ in_data.fast_pairs[i][0] ][j] = 0;
        }
    }
    
    
    /*
    ============= Need to base this on averaging =============
    */
    
    // Choose a reaction to fire, state to fire from, and time step
    
    
    r_timestep = ((double) rand() / (RAND_MAX));
    dt = log(1 / r_timestep) / props[2];
    
    rxn_to_fire_ind = 2;
    
    int AB_total = N[0] + N[1];
    N[0] = (int) N_micro_avg[0];
    N[1] = AB_total - N[0];

    /*
    ============= Sensitivity analysis part =============
    */
    
    // Microscale sensitivities
    micro_scale_sens[0][0] = - N_micro_avg[0] / (in_data.rate_const[0] + in_data.rate_const[1]);
    micro_scale_sens[0][1] = N_micro_avg[1] / (in_data.rate_const[0] + in_data.rate_const[1]);
    micro_scale_sens[0][2] = 0;
    micro_scale_sens[1][0] = N_micro_avg[0] / (in_data.rate_const[0] + in_data.rate_const[1]);
    micro_scale_sens[1][1] = - N_micro_avg[1] / (in_data.rate_const[0] + in_data.rate_const[1]);
    micro_scale_sens[1][2] = 0;
    micro_scale_sens[2][0] = 0;
    micro_scale_sens[2][1] = 0;
    micro_scale_sens[2][2] = 0;
    
    // Derivative of propensities - indirect
    prop_ders_indirect[0][0] = 0;
    prop_ders_indirect[0][1] = 0;
    prop_ders_indirect[0][2] = 0;
    prop_ders_indirect[1][0] = 0;
    prop_ders_indirect[1][1] = 0;
    prop_ders_indirect[1][2] = 0;
    prop_ders_indirect[2][0] = in_data.rate_const[2] * N_micro_avg[1] / (in_data.rate_const[0] + in_data.rate_const[1]);
    prop_ders_indirect[2][1] = - in_data.rate_const[2] * N_micro_avg[1] / (in_data.rate_const[0] + in_data.rate_const[1]);
    prop_ders_indirect[2][2] = 0;
    
    
    // Derivative of propensities - total (sum of direct and indirect)
    for (int i = 0; i < in_data.n_rxns; i++){
        for (int j = 0; j < in_data.n_params; j++){
            prop_ders[i][j] = prop_ders_direct[i][j] + prop_ders_indirect[i][j];
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