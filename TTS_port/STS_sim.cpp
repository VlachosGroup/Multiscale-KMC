# include <iostream>    // using IO functions
# include <string>      // using string
# include <vector>
# include <fstream>
# include <random>
# include <math.h>

# include <typeinfo>
# include <sstream>
# include <stdio.h>      /* printf, fgets */
# include <stdlib.h>     /* atoi */
# include "myHeader.h"
using namespace std;


string STS_traj :: species_out_flname = "specs.out";
string STS_traj :: traj_deriv_out_flname = "SA.out";

STS_traj :: STS_traj() : in_data(){}

void STS_traj :: simulate(int rand_seed){

cout << "Trajectory running KMC simulation." << endl;

/*
============================ Initialize variables ============================
*/



    
//srand(rand_seed);      // Set the random seed
//
//
//ind_rec_spec = 0;
//ind_rec_rxns = 0;
//
//t = 0;
//t_prev = 0;
//ind_rec = 0;
//
//// Initialize species populations
//for(int k =0; k < n_specs; k++){
//    N[k] = N_0[k];
//}
//
//// Initialize trajectory derivatives
//for(int k =0; k < n_params; k++){
//    W[k] = 0;
//}
//
//// Open files to record data for the trajectory
//ofstream writer_spec;
//ofstream writer_SA;
//
//if(write_traj_files){
//    
//    // Create a folder for all trajectory data files
//    // Name the file specific to this trajectory. Put random seed in the file name.
//    
//    writer_spec.open(STS_traj.species_out_flname);
//    if(! writer_spec){  
//        cout << "Error opening file" << endl;
//        return -1;
//    }
//    
//    // Write header    
//    writer_spec << "Time \t";
//    for(int j = 0; j < n_specs; j++){
//        writer_spec << spec_names[j] << "\t";
//    }
//    writer_spec << endl;
//    
//    writer_SA.open(STS_traj.species_out_flname);
//    if(! writer_SA){  
//        cout << "Error opening file" << endl;
//        return -1;
//    }
//    
//    // Write header
//    writer_SA << "Time \t";
//    for(int i = 0; i < n_params; i++){
//        writer_SA << param_names[i] << "\t";
//    }
//    writer_SA << endl; 
//}        
//
//// Start KMC loop
//while(t < t_final){
//    
//    /*
//    ============================ Compute quantities for current step ============================
//    */
//    
//    // Generate random numbers
//    r_rxn_choose = ((double) rand() / (RAND_MAX));
//    r_timestep = ((double) rand() / (RAND_MAX));
//    
//    // Compute reaction propensities
//    for(int i = 0; i < n_rxns; i++){
//        
//        // Use mass action kinetics equation, a = k * [A]^ma * [B]^mb * ...
//        props[i] = rate_const[i];
//        for(int j = 0; j < n_specs; j++){
//            if(stoich_mat[i][j] < 0){//if this species is a reactant, use it to compute the rate
//                props[i] = props[i] * pow (N[j], -stoich_mat[i][j]);
//            }
//        }
//        
//        // Fill in derivatives for each parameter
//        for(int k = 0; k < n_params; k++){
//            if(i==k){
//                
//                prop_ders[i][k] = 1.0;
//                
//                for(int j = 0; j < n_specs; j++){
//                    if(stoich_mat[i][j] < 0){//if this species is a reactant, use it to compute the rate
//                        prop_ders[i][k] = prop_ders[i][k] * pow (N[j], -stoich_mat[i][j]);
//                    }
//                }
//                
//            }else{
//                prop_ders[i][k] = 0;
//            }
//        }
//    }
//    
//    // Add derivatives for each parameter
//    for(int k=0; k < n_params; k++){
//        prop_ders_sum[k] = 0;
//        for(int i=0; i < n_rxns; i++){
//            prop_ders_sum[k] += prop_ders[i][k];
//        }
//    }
//    
//    // Sum the propensities
//    asum = 0;
//    for(int i = 0; i < n_rxns; i++){
//        asum += props[i];
//    }
//    
//    // If all propensities are 0, then exit the while loop
//    if(asum == 0){
//        break;
//    }
//    
//    // Choose which reaction to fire
//    for(int i = 0; i < n_rxns; i++){
//        prop_cum[i] = 0;
//        for(int j = 0; j <= i; j++){
//            prop_cum[i] += props[j] / asum;
//        }
//    }
//    
//    rxn_to_fire_ind = 0;
//    while(prop_cum[rxn_to_fire_ind] < r_rxn_choose){
//        rxn_to_fire_ind += 1;
//    }
//    
//    // Compute time step
//    dt = log(1 / r_timestep) / asum;
//    
//    // Record the current state as long as time >= t_sample
//    while(t >= t_rec[ind_rec]){
//        
//        record_data()
//        
//    }
//    
//    /*
//    ============ Fire reaction and update system ==============
//    */
//    
//    // Update species populations
//    for(int j = 0; j < n_specs; j++){
//        N[j] += stoich_mat[rxn_to_fire_ind][j];
//    }
//    
//    // Update clock
//    t_prev = t;
//    t = t + dt;
//    
//    // Update trajectory derivatives
//    for(int k=0; k < n_params; k++){
//        
//        W[k] += prop_ders[rxn_to_fire_ind][k] / props[rxn_to_fire_ind];         // contribution from reaction that fires            
//        W[k] -= prop_ders_sum[k] * dt;       // contribution from time step
//    }
//
//}
//
//// Fill in recording times that were missed
//while(ind_rec <= N_record){
//    
//	record_data();
//        
//    }
//
//// Close output files
//if(write_traj_files){
//    writer_spec.close();
//    writer_SA.close();
//}

}


/*
========== Record data in the KMC simulation ==========
*/
void STS_traj :: record_data(){

//    t_trunc = t_rec[ind_rec] - t_prev;
//        
//    // Record data in output files
//    if(write_traj_files){
//        
//        // Record the species numbers
//        writer_spec << t_rec[ind_rec] << "\t";
//        for(int j = 0; j < n_specs; j++){
//            writer_spec << N[j] << "\t";
//        }
//        writer_spec << endl;
//    
//        // Record trajectory derivatives
//        writer_SA << t_rec[ind_rec] << "\t";
//        for(int k = 0; k < n_params; k++){
//            writer_SA << W[k] - prop_ders_sum[k] * t_trunc << "\t";
//        }
//        writer_SA << endl;
//    }
//    
//    // Record data in variables
//    
//    // Record species populations
//    for(int j = 0; j < n_specs; j++){
//         spec_profiles_pp[ind_rec_spec] = N[j];
//         ind_rec_spec += 1;
//        }
//    
//    // Record trajectory derivatives
//    for(int k = 0; k < n_params; k++){
//        traj_derivs_pp[ind_rec_rxns] = W[k] - prop_ders_sum[k] * t_trunc;
//        ind_rec_rxns += 1;
//        }
//    
//    
//    ind_rec += 1;
	
}