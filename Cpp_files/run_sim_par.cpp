// Parallel (MPI) version of running many KMC simulations

# include <mpi.h>
# include <iostream>    // using IO functions
# include <string>      // using string
# include <vector>
# include <fstream>
# include <random>
# include <math.h>
# include "myHeader.h"
using namespace std;

void Traj_stats :: run_simulations(){       // add some if statements which will change bahavior depending on if it is STS or TTS

    /*
    ============================ Run simulations to gather data ============================
    */
    
    initialize_stats();

    // Set up MPI variables
    
    int id;
    int ierr;
    int p;
    int Npp;
    //int argc;
    //char *argv[];
    
    ierr = MPI_Init ( NULL, NULL );                 //  Initialize MPI.
    ierr = MPI_Comm_size ( MPI_COMM_WORLD, &p );      //  Get the number of processes.  
    ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id );     //  Get the individual process ID.
    
    
    
    Npp = in_data.N_traj / p + 1;
    in_data.N_traj = Npp * p;       // round up to the nearest multiple of the number of processors
    
    if(id==0){
        cout << in_data.N_traj << " replicate trajectories will be used." << endl;
    }
    
    
    // Make a data array for each processor to be combined later
    
    int size1 = in_data.N_record * in_data.n_specs;
    int size2 = in_data.N_record * in_data.n_params;
    int size3 = in_data.N_record * in_data.n_specs * in_data.n_params;
    
    double spec_profiles_averages_arr[size1];
    double traj_deriv_avgs_arr[size2];
    double sensitivities_arr[size3];
    
    for(int i = 0; i < size1; i++){
        spec_profiles_averages_arr[i] = 0;
    }
    
    for(int i = 0; i < size2; i++){
        traj_deriv_avgs_arr[i] = 0;
    }
    
    for(int i = 0; i < size3; i++){
        sensitivities_arr[i] = 0;
    }
    
    int ind1;
    int ind2;
    int ind3;
    
    /*
    ============ run KMC simulations to gather data ==============
    */       

    for(int traj_ind = 0; traj_ind < Npp; traj_ind++){
        
        // Create and run a KMC simulation
        KMC_traj run;           // Change this like for TTS - make TTS_traj object instead
        run.in_data = in_data;      // Copy input file data to the trajectory object
        
        if(! in_data.two_time_scale){
            run.simulate_STS(12345  + id * Npp + traj_ind);
        }else{
            run.simulate_TTS(12345  + id * Npp + traj_ind);
        }
        
        
        // Add to statistical running counts
        
        ind1 = 0;
        ind2 = 0;
        ind3 = 0;
        
        for (int i = 0; i < in_data.N_record; ++i){
            
            for(int j = 0; j < in_data.n_specs; j++){
                spec_profiles_averages_arr[ind1] += run.spec_profile[i][j];
                ind1++;
                
                for(int k = 0; k < in_data.n_params; k++){
                    sensitivities_arr[ind3] += run.spec_profile[i][j] * run.traj_deriv_profile[i][k];    // Change this like for TTS - add microscale contribution
                    
                    ind3++;
                }
            }
            
            for(int j = 0; j < in_data.n_params; j++){
                traj_deriv_avgs_arr[ind2] += run.traj_deriv_profile[i][j];
                ind2++;
            }
        }
        
        
    }

    /*
    ============ Collect data from the replicate simulations ==============
    */
    
    // For final analysis, gather data into the big arrays
    double spec_profiles_averages_big_arr[p * size1];
    double traj_deriv_avgs_big_arr[p * size2];
    double sensitivities_big_arr[p * size3];
    
    MPI_Gather(spec_profiles_averages_arr, size1, MPI_DOUBLE, spec_profiles_averages_big_arr, size1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(traj_deriv_avgs_arr, size2, MPI_DOUBLE, traj_deriv_avgs_big_arr, size2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(sensitivities_arr, size3, MPI_DOUBLE, sensitivities_big_arr, size3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    
    // Reshape the data into multidimensional arrays on processor 0
    if(id==0){
        
        ind1 = 0;
        ind2 = 0;
        ind3 = 0;
        
        //for 
        
        // Normalize statistics by the number of trajectories
        
        for(int proc_id = 0; proc_id < p; proc_id++){
        
            for (int i = 0; i < in_data.N_record; ++i){
                    
                for(int j = 0; j < in_data.n_specs; j++){
                    
                    spec_profiles_averages[i][j] += spec_profiles_averages_big_arr[ind1];
                    ind1++;
                    
                    for(int k = 0; k < in_data.n_params; k++){
                        
                        sensitivities[i][j][k] += sensitivities_big_arr[ind3];
                        ind3++;
                        
                    }
                }
                
                for(int j = 0; j < in_data.n_params; j++){
                        traj_deriv_avgs[i][j] += traj_deriv_avgs_big_arr[ind2];
                        ind2++;
                }
                    
                    
            }
        }
        
        
        finalize_stats();
        write_spec_avg_output();                // Write output files   
        write_sensitivity_output();
        cout << "Simulation completed successfully" << endl;
    }
    
    MPI_Finalize ( );     //  Terminate MPI.
    
}