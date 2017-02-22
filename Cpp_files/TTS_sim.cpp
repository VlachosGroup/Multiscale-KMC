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




}


/*
============= Microscale KMC simulation =============
*/

void KMC_traj_TTS :: simulate_micro(){      // Implement with analytical solution for A->B for now

    //double r_rxn_choose;                    // random number between 0 and 1 used to choose which reaction to fire
    //double r_timestep;                      // random mumber between 0 and 1 used to choose the time step
    //
    ////int micro_step_limit = 10000;
    ////for(int micro_step = 0, micro_step < micro_step_limit, micro_step++){
    ////    // Execute metropolis simulation
    ////}
    //
    //int AB_total = N[0] + N[1];
    //double K = in_data.rate_const[0] / in_data.rate_const[1];
    //
    //// Average species
    //N_micro_avg[0] = AB_total * 1 / (K + 1);
    //N_micro_avg[1] = AB_total * K / (K + 1);
    //N_micro_avg[2] = N[2];
    //
    //// Average propensities
    //
    //
    //// Derivative of propensities
    
    
}


/*
========= Record the current state
*/
void KMC_traj_TTS :: record_stats_TTS(){
    
    //t_trunc = in_data.t_rec[ind_rec] - t_prev;
    //    
    //// Record data in output files
    //if(in_data.write_traj_files){
    //    
    //    // Record the species numbers
    //    writer_spec << in_data.t_rec[ind_rec] << "\t";
    //    for(int j = 0; j < in_data.n_specs; j++){
    //        writer_spec << N_micro_avg[j] << "\t";
    //    }
    //    writer_spec << endl;
    //
    //    // Record trajectory derivatives
    //    writer_SA << in_data.t_rec[ind_rec] << "\t";
    //    for(int k = 0; k < in_data.n_params; k++){
    //        writer_SA << W[k] - prop_ders_sum[k] * t_trunc << "\t";
    //    }
    //    writer_SA << endl;
    //}
    //
    //// Record data in variables
    //
    //// Record species populations
    //for(int j = 0; j < in_data.n_specs; j++){
    //    spec_profile[ind_rec][j] = N_micro_avg[j];
    //    }
    //
    //// Record trajectory derivatives
    //for(int k = 0; k < in_data.n_params; k++){
    //    traj_deriv_profile[ind_rec][k] = W[k] - prop_ders_sum[k] * t_trunc;
    //    }
    //
    //ind_rec += 1;
    
}