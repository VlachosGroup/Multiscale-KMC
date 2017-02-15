# include <iostream>    // using IO functions
# include <string>      // using string
# include <vector>
# include <fstream>
# include <random>
# include <math.h>
# include <sstream>
# include "myHeader.h"
using namespace std;

// Static variables - names of output files
string STS_traj :: species_out_flname = "specs";
string STS_traj :: traj_deriv_out_flname = "SA";

STS_traj :: STS_traj() : in_data(){}            // empty constructor

void STS_traj :: simulate(int rand_seed){       // Execute simulation


/*
============================ Initialize variables ============================
*/

double props[in_data.n_rxns];
double prop_cum[in_data.n_rxns];
double asum;

// Used for sensitivity analysis
double W[in_data.n_params];                       // trajectory derivatives
double prop_ders[in_data.n_rxns][in_data.n_params];     // derivative of each propensity with respect to each parameter
double prop_ders_sum[in_data.n_params];           // derivative of the sum of all propensities with respect to each parameter
double t_trunc;                         // time since the previous KMC step

double t = 0;                               // KMC clock
double t_prev = 0;                          // KMC time of previous step
double dt;                              // time step
double r_rxn_choose;                    // random number between 0 and 1 used to choose which reaction to fire
double r_timestep;                      // random mumber between 0 and 1 used to choose the time step
int rxn_to_fire_ind;                    // index of the reaction chosen to fire
int ind_rec = 0;                            // time point

// Vector for recording species profiles
spec_profile.resize(in_data.N_record);
for (int i = 0; i < in_data.N_record; i++){
    spec_profile[i].resize(in_data.n_specs);}

// Vector for recording trajectory derivatives
traj_deriv_profile.resize(in_data.N_record);
for (int i = 0; i < in_data.N_record; i++){
    traj_deriv_profile[i].resize(in_data.n_params);}
    
srand(rand_seed);      // Set the random seed


// Initialize species populations
int N[in_data.n_specs];

for(int k =0; k < in_data.n_specs; k++){
    N[k] = in_data.N_0[k];
}

// Initialize trajectory derivatives
for(int k =0; k < in_data.n_params; k++){
    W[k] = 0;
}

/*
============== Open files to record data for the trajectory ==============
*/
ofstream writer_spec;
ofstream writer_SA;

if(in_data.write_traj_files){
    
    // Create a folder for all trajectory data files
    // Name the file specific to this trajectory. Put random seed in the file name.
    
    ostringstream fname1;
    fname1 << STS_traj::species_out_flname << "_" << rand_seed << ".out";
    
    writer_spec.open(fname1.str());
    if(! writer_spec){  
        cout << "Error opening file" << endl;
    }
    
    // Write header    
    writer_spec << "Time \t";
    for(int j = 0; j < in_data.n_specs; j++){
        writer_spec << in_data.spec_names[j] << "\t";
    }
    writer_spec << endl;
    
    ostringstream fname2;
    fname2 << STS_traj::traj_deriv_out_flname << "_" << rand_seed << ".out";
    
    writer_SA.open(fname2.str());
    if(! writer_SA){  
        cout << "Error opening file" << endl;
    }
    
    // Write header
    writer_SA << "Time \t";
    for(int i = 0; i < in_data.n_params; i++){
        writer_SA << in_data.param_names[i] << "\t";
    }
    writer_SA << endl; 
}        

// Start KMC loop
while(t < in_data.t_final){
    
    /*
    ============================ Compute quantities for current step ============================
    */
    
    // Generate random numbers
    r_rxn_choose = ((double) rand() / (RAND_MAX));
    r_timestep = ((double) rand() / (RAND_MAX));
    
    // Compute reaction propensities
    for(int i = 0; i < in_data.n_rxns; i++){
        
        // Use mass action kinetics equation, a = k * [A]^ma * [B]^mb * ...
        props[i] = in_data.rate_const[i];
        for(int j = 0; j < in_data.n_specs; j++){
            if(in_data.stoich_mat[i][j] < 0){//if this species is a reactant, use it to compute the rate
                props[i] = props[i] * pow (N[j], -in_data.stoich_mat[i][j]);
            }
        }
        
        // Fill in derivatives for each parameter
        for(int k = 0; k < in_data.n_params; k++){
            if(i==k){
                
                prop_ders[i][k] = 1.0;
                
                for(int j = 0; j < in_data.n_specs; j++){
                    if(in_data.stoich_mat[i][j] < 0){//if this species is a reactant, use it to compute the rate
                        prop_ders[i][k] = prop_ders[i][k] * pow (N[j], -in_data.stoich_mat[i][j]);
                    }
                }
                
            }else{
                prop_ders[i][k] = 0;
            }
        }
    }
    
    // Add derivatives for each parameter
    for(int k=0; k < in_data.n_params; k++){
        prop_ders_sum[k] = 0;
        for(int i=0; i < in_data.n_rxns; i++){
            prop_ders_sum[k] += prop_ders[i][k];
        }
    }
    
    // Sum the propensities
    asum = 0;
    for(int i = 0; i < in_data.n_rxns; i++){
        asum += props[i];
    }
    
    // If all propensities are 0, then exit the while loop
    if(asum == 0){
        break;
    }
    
    // Choose which reaction to fire
    for(int i = 0; i < in_data.n_rxns; i++){
        prop_cum[i] = 0;
        for(int j = 0; j <= i; j++){
            prop_cum[i] += props[j] / asum;
        }
    }
    
    rxn_to_fire_ind = 0;
    while(prop_cum[rxn_to_fire_ind] < r_rxn_choose){
        rxn_to_fire_ind += 1;
    }
    
    // Compute time step
    dt = log(1 / r_timestep) / asum;
    
    /*
    ========= Record the current state as long as time >= t_sample
    */
    while(t >= in_data.t_rec[ind_rec]){
        
        t_trunc = in_data.t_rec[ind_rec] - t_prev;
        
        // Record data in output files
        if(in_data.write_traj_files){
            
            // Record the species numbers
            writer_spec << in_data.t_rec[ind_rec] << "\t";
            for(int j = 0; j < in_data.n_specs; j++){
                writer_spec << N[j] << "\t";
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
            spec_profile[ind_rec][j] = (double) N[j];
            }
        
        // Record trajectory derivatives
        for(int k = 0; k < in_data.n_params; k++){
            traj_deriv_profile[ind_rec][k] = W[k] - prop_ders_sum[k] * t_trunc;
            }
        
        ind_rec += 1;
        
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
    
    t_trunc = in_data.t_rec[ind_rec] - t_prev;
        
    // Record data in output files
    if(in_data.write_traj_files){
        
        // Record the species numbers
        writer_spec << in_data.t_rec[ind_rec] << "\t";
        for(int j = 0; j < in_data.n_specs; j++){
            writer_spec << N[j] << "\t";
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
        spec_profile[ind_rec][j] = N[j];
        }
    
    // Record trajectory derivatives
    for(int k = 0; k < in_data.n_params; k++){
       traj_deriv_profile[ind_rec][k] = W[k] - prop_ders_sum[k] * t_trunc;
        }
    
    ind_rec += 1;
        
    }

// Close output files
if(in_data.write_traj_files){
    writer_spec.close();
    writer_SA.close();
}


}

