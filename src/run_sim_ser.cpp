// Serial version of running many KMC simulations

# include <vector>
# include <random>
# include <math.h>
# include "myHeader.h"
using namespace std;

/*
============================ Run simulations to gather data ============================
*/

void Traj_stats :: run_simulations(){       // add some if statements which will change bahavior depending on if it is STS or TTS

    initialize_stats();   

    for(int traj_ind = 0; traj_ind < in_data.N_traj; traj_ind++){
        
        // Create and run a KMC simulation
        cout << traj_ind << " / " << in_data.N_traj << endl;
        
        if(! in_data.two_time_scale){       // Single time scale
            
            KMC_traj run;           // Change this like for TTS - make TTS_traj object instead
            run.in_data = in_data;      // Copy input file data to the trajectory object
            run.simulate(12345  + traj_ind);
            
            // Add to statistical running counts
            for (int i = 0; i < in_data.N_record; ++i){
                
                for(int j = 0; j < in_data.n_specs; j++){
                    spec_profiles_averages[i][j] += run.spec_profile[i][j];
                    
                    for(int k = 0; k < in_data.n_params; k++){
                        sensitivities[i][j][k] += run.spec_profile[i][j] * run.traj_deriv_profile[i][k];
                    }
                }
                
                for(int j = 0; j < in_data.n_params; j++){
                    traj_deriv_avgs[i][j] += run.traj_deriv_profile[i][j];
                }
            }
            
        }else{                          // Two time scale
            
            KMC_traj_TTS run;           // Change this like for TTS - make TTS_traj object instead
            run.in_data = in_data;      // Copy input file data to the trajectory object
            run.simulate_TTS(12345  + traj_ind);
            
            // Add to statistical running counts
            for (int i = 0; i < in_data.N_record; ++i){
                
                for(int j = 0; j < in_data.n_specs; j++){
                    spec_profiles_averages[i][j] += run.spec_profile[i][j];
                    
                    for(int k = 0; k < in_data.n_params; k++){
                        sensitivities[i][j][k] += run.spec_profile[i][j] * run.traj_deriv_profile[i][k] + run.micro_scale_sens_profile[i][j][k];    // Add microscale contribution
                    }
                }
                
                for(int j = 0; j < in_data.n_params; j++){
                    traj_deriv_avgs[i][j] += run.traj_deriv_profile[i][j];
                }
            }
            
        }
        
        
    }
    

    finalize_stats();

    // Write output files
    write_spec_avg_output();        
    write_sensitivity_output();
   
    cout << "Simulation complete." << endl;
}