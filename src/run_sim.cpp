// Serial version of running many KMC simulations

# include <vector>
# include <random>
# include <math.h>
# include "myHeader.h"
using namespace std;

/*
============================ Run simulations to gather data ============================
*/

#ifndef MPI

void Traj_stats :: run_simulations(){

    initialize_stats();

    for(int traj_ind = 0; traj_ind < in_data.N_traj; traj_ind++){

        // Create and run a KMC simulation

        //cout << traj_ind + 1 << " / " << in_data.N_traj << endl;              // Print the trajectory number

        KMC_traj* run = NULL;  // initalize the pointer

        if(in_data.two_time_scale){     // Two time scale
            run = new KMC_traj_TTS;
        }else{                          // Single time scale
            run = new KMC_traj;
        }

        run->in_data = in_data;             // Copy input file data to the trajectory object
        run->simulate(in_data.rand_seed + traj_ind);

        // Add to statistical running counts
        for (int i = 0; i < in_data.N_record; ++i){

            for(int j = 0; j < in_data.n_specs; j++){
                spec_profiles_averages[i][j] += run->spec_profile[i][j];

                for(int k = 0; k < in_data.n_params; k++){  // For TTS, add extra contribution from microscale averaging
                    sensitivities[i][j][k] += run->spec_profile[i][j] * run->traj_deriv_profile[i][k] + run->get_micro_scale_sens_profile(i, j, k);
                        if ( not std::isfinite( sensitivities[i][j][k] ) ){
                            cout << traj_ind << " has NaNs" << endl;
                            cout << in_data.rand_seed + traj_ind << " is the random seed." << endl;
                        }
                }
            }

            for(int j = 0; j < in_data.n_params; j++){
                traj_deriv_avgs[i][j] += run->traj_deriv_profile[i][j];
            }
        }

    }


    finalize_stats();

    // Write output files
    write_spec_avg_output();
    write_sensitivity_output();

    cout << "Simulation complete." << endl;
}

#endif
