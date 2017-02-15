// Serial version of running many KMC simulations

# include <iostream>    // using IO functions
# include <string>      // using string
# include <vector>
# include <fstream>
# include <random>
# include <math.h>
# include "myHeader.h"
using namespace std;

/*
============================ Run simulations to gather data ============================
*/

void Traj_stats_STS :: run_simulations(){

    initialize_stats();   

    for(int traj_ind = 0; traj_ind < in_data.N_traj; traj_ind++){
        
        // Create and run a KMC simulation
        STS_traj run;           // Change this like for TTS - make TTS_traj object instead
        run.in_data = in_data;      // Copy input file data to the trajectory object
        run.simulate(12345  + traj_ind);
        
        // Add to statistical running counts
        for (int i = 0; i < in_data.N_record; ++i){
            
            for(int j = 0; j < in_data.n_specs; j++){
                spec_profiles_averages[i][j] += run.spec_profile[i][j];
                
                for(int k = 0; k < in_data.n_params; k++){
                    sensitivities[i][j][k] += run.spec_profile[i][j] * run.traj_deriv_profile[i][k];    // Change this like for TTS - add microscale contribution
                }
            }
            
            for(int j = 0; j < in_data.n_params; j++){
                traj_deriv_avgs[i][j] += run.traj_deriv_profile[i][j];
            }
        }
    }
    

    finalize_stats();

    // Write output files
    write_spec_avg_output();        
    write_sensitivity_output();
   
    cout << "Simulation complete." << endl;
}


// void Traj_stats_TTS :: run_simulations(){}       // Need to set up this method to be overwritten by child class