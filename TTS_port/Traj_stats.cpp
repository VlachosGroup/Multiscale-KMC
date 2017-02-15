# include <iostream>    // using IO functions
# include <string>      // using string
# include <vector>
# include <fstream>
# include <random>
# include <math.h>
# include "myHeader.h"
using namespace std;

string Traj_stats_STS :: species_avgs_out_flname = "species_avgs.out";
string Traj_stats_STS :: SA_out_flname = "sensitivities.out";

Traj_stats_STS :: Traj_stats_STS() : in_data(){}      // Empty class constructor


void Traj_stats_STS :: initialize_stats(){
    
    // Resize vectors to be able to hold the necessary data
    
    spec_profiles_averages.resize(in_data.N_record);
    for (int i = 0; i < in_data.N_record; i++){
        spec_profiles_averages[i].resize(in_data.n_specs);}

    traj_deriv_avgs.resize(in_data.N_record);
    for (int i = 0; i < in_data.N_record; i++){
        traj_deriv_avgs[i].resize(in_data.n_params);}
        
    sensitivities.resize(in_data.N_record);
    for (int i = 0; i < in_data.N_record; i++){
        
        sensitivities[i].resize(in_data.n_specs);
        for(int j = 0; j < in_data.n_specs; j++){
            
            sensitivities[i][j].resize(in_data.n_params);
        }
    }


    // Intialize statistics
    
    for (int i = 0; i < in_data.N_record; ++i){
        for(int j = 0; j < in_data.n_specs; j++){
            spec_profiles_averages[i][j] = 0;
        }
    }

    for (int i = 0; i < in_data.N_record; ++i){
        for(int j = 0; j < in_data.n_params; j++){
            traj_deriv_avgs[i][j] = 0;
        }
    }
    
    for (int i = 0; i < in_data.N_record; ++i){
        for(int j = 0; j < in_data.n_specs; j ++){
            for(int k = 0; k < in_data.n_params; k++){
                sensitivities[i][j][k] = 0;
            }
        }
    }
    
}


void Traj_stats_STS :: finalize_stats(){

    /*
    ============ Finalize statistical analysis ==============
    */

    // Normalize statistics by the number of trajectories
    for (int i = 0; i < in_data.N_record; ++i){
            
        for(int j = 0; j < in_data.n_specs; j++){
            spec_profiles_averages[i][j] = spec_profiles_averages[i][j] / in_data.N_traj;
            
            for(int k = 0; k < in_data.n_params; k++){
                sensitivities[i][j][k] = sensitivities[i][j][k] / in_data.N_traj;
            }
        }
        
        for(int j = 0; j < in_data.n_params; j++){
            traj_deriv_avgs[i][j] = traj_deriv_avgs[i][j] / in_data.N_traj;
        }
    }
    
    
    // Subtract 2nd term of the covariance, renormalize by N_traj-1
    for (int i = 0; i < in_data.N_record; ++i){
            
        for(int j = 0; j < in_data.n_specs; j++){
            
            for(int k = 0; k < in_data.n_params; k++){
                sensitivities[i][j][k] -= spec_profiles_averages[i][j] * traj_deriv_avgs[i][k];
                sensitivities[i][j][k] = sensitivities[i][j][k] * in_data.N_traj / (in_data.N_traj - 1);
            }
        }
    }
   
}

/* 
=========== Print species population averages into an output file ===========
*/
void Traj_stats_STS :: write_spec_avg_output(){

    ofstream writer_spec_avgs;
    
    writer_spec_avgs.open(Traj_stats_STS :: species_avgs_out_flname);
    if(! writer_spec_avgs){  
        cout << "Error opening file" << endl;
    }
    
    // Write header
    writer_spec_avgs << "Species mean populations versus time" << endl;
    writer_spec_avgs << endl;
    writer_spec_avgs << "Time \t";
    for(int j = 0; j < in_data.n_specs; j++){
        writer_spec_avgs << in_data.spec_names[j] << "\t";
    }
    writer_spec_avgs << endl;
    
    for(int time_ind = 0; time_ind < in_data.N_record; time_ind ++){
        
        writer_spec_avgs << in_data.t_rec[time_ind] << "\t";
        
        for(int spec_ind = 0; spec_ind < in_data.n_specs; spec_ind++){
            writer_spec_avgs << spec_profiles_averages[time_ind][spec_ind] << "\t"; // print the mean population of this species at this time
        }
        
        writer_spec_avgs << endl;
    }
    
    writer_spec_avgs.close();

}
        

/*
===========Print sensitivities into an output file ===========
*/
void Traj_stats_STS :: write_sensitivity_output(){
    
    ofstream writer_sensitivities;

    writer_sensitivities.open(Traj_stats_STS :: SA_out_flname);
    if(! writer_sensitivities){  
        cout << "Error opening file" << endl;
    }
    
    // Write header
    writer_sensitivities << "Sensitivities for each species and rate constant versus time" << endl << endl;
    
    // Loop over species
    for(int i=0; i < in_data.n_specs; i++){
        
        writer_sensitivities << in_data.spec_names[i] << endl << endl;
        
        // Write header for this species
        writer_sensitivities << "Time \t";
        for(int k = 0; k < in_data.n_params; k++){
            writer_sensitivities << in_data.param_names[k] << "\t";
        }
        writer_sensitivities << endl;
        
        // Loop over parameters to print out 
        for(int time_ind = 0; time_ind < in_data.N_record; time_ind ++){
            
            writer_sensitivities << in_data.t_rec[time_ind] << "\t";
            
            for(int param_ind = 0; param_ind < in_data.n_params; param_ind++){
                writer_sensitivities << sensitivities[time_ind][i][param_ind] << "\t";
            }
            
            writer_sensitivities << endl;
        }
        
        // Put a barrier in between species data
        if(i < in_data.n_specs - 1){
            writer_sensitivities << endl;
            writer_sensitivities << "=======================================" << endl;
            writer_sensitivities << endl;
        }
        
    }
    
    writer_sensitivities.close();

}
    

// Main function
int main() {
	
    cout.precision(6);
	file_reader fr("network.in");               // Read input file
        
    if(! fr.two_time_scale){                    // single time scale
            
        Traj_stats_STS STS_sim; 
        STS_sim.in_data = fr;                   // assign input
        STS_sim.run_simulations();              // run KMC simulations
        
    }else{                                      // two time scale version
        
        Traj_stats_TTS TTS_sim;
        TTS_sim.in_data = fr;
        TTS_sim.run_simulations();
        
    }
    
    

	return 0;
}