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

    //double r_rxn_choose;                    // random number between 0 and 1 used to choose which reaction to fire
    //double r_accept;                      // random mumber between 0 and 1 used to accept or reject the reaction
    
    //int micro_step_limit = 10000;
    //for(int micro_step = 0, micro_step < micro_step_limit, micro_step++){
    //    // Execute metropolis simulation
    //}

    
    double prop_ders_direct[in_data.n_rxns][in_data.n_params];
    double prop_ders_indirect[in_data.n_rxns][in_data.n_params];
    
    int AB_total = N[0] + N[1];
    double K = in_data.rate_const[0] / in_data.rate_const[1];
    
    // Average species
    N_micro_avg[0] = AB_total * 1 / (K + 1);
    N_micro_avg[1] = AB_total * K / (K + 1);
    N_micro_avg[2] = N[2];
    
    // Average propensities
    
    props[0] = 0;
    props[1] = 0;
    props[2] = in_data.rate_const[2] * N_micro_avg[1];
    
    double r_timestep;
    r_timestep = ((double) rand() / (RAND_MAX));
    dt = log(1 / r_timestep) / props[2];
    
    rxn_to_fire_ind = 2;
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
    
    // Derivative of propensities - direct
    prop_ders_indirect[0][0] = 0;
    prop_ders_indirect[0][1] = 0;
    prop_ders_indirect[0][2] = 0;
    prop_ders_indirect[1][0] = 0;
    prop_ders_indirect[1][1] = 0;
    prop_ders_indirect[1][2] = 0;
    prop_ders_indirect[2][0] = in_data.rate_const[2] * N_micro_avg[1] / (in_data.rate_const[0] + in_data.rate_const[1]);
    prop_ders_indirect[2][1] = - in_data.rate_const[2] * N_micro_avg[1] / (in_data.rate_const[0] + in_data.rate_const[1]);
    prop_ders_indirect[2][2] = 0;
    
    
    // Derivative of propensities - indirect
    prop_ders_direct[0][0] = 0;
    prop_ders_direct[0][1] = 0;
    prop_ders_direct[0][2] = 0;
    prop_ders_direct[1][0] = 0;
    prop_ders_direct[1][1] = 0;
    prop_ders_direct[1][2] = 0;
    prop_ders_direct[2][0] = 0;
    prop_ders_direct[2][1] = 0;
    prop_ders_direct[2][2] = N_micro_avg[1];
    
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