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

void Traj_stats_STS :: run_simulations(){

    /*
    ============================ Run simulations to gather data ============================
    */
    
    

    // Set up MPI variables
    
    int id;
    int ierr;
    int p;
    int Npp;
    
    ierr = MPI_Init ( &argc, &argv );                 //  Initialize MPI.
    ierr = MPI_Comm_size ( MPI_COMM_WORLD, &p );      //  Get the number of processes.  
    ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id );     //  Get the individual process ID.
    
    
    
    Npp = N_traj / p + 1;
    N_traj = Npp * p;       // round up to the nearest multiple of the number of processors
    
    if(id==0){
        cout << N_traj << " replicate trajectories will be used." << endl;
    }
    
    
    int seed_list[N_traj];
    for(int i = 0; i < N_traj; i++){
        seed_list[i] = i+1;
    }
    
    
    // For each processor
    int spec_ints_pp = Npp * (N_record+1) * n_specs;
    int traj_ints_pp = Npp * (N_record+1) * n_params;
    int spec_profiles_pp[spec_ints_pp];
    double traj_derivs_pp[traj_ints_pp];
    // Fill these up with data as the simulation progresses
    
    // Scatter the seed list
    int sub_seeds[Npp];
    MPI_Scatter(seed_list, Npp, MPI_INT, sub_seeds,
            Npp, MPI_INT, 0, MPI_COMM_WORLD);
    


    double spec_avgs[in_data.N_record][in_data.n_specs];
    double W_avgs[in_data.N_record][in_data.n_params];
    double prod_avgs[in_data.N_record][in_data.n_specs][in_data.n_params];

    /*
    ============ run KMC simulations to gather data ==============
    */       

    STS_traj run;
    run.in_data = in_data;      // Copy input file data to the trajectory object
    run.simulate(12345);

    /*
    ============ Collect data from the replicate simulations ==============
    */
    
    // For final analysis, gather data into the big arrays
    int spec_profiles[N_traj * (N_record+1) * n_specs];
    double traj_derivs[N_traj * (N_record+1) * n_params];
    
    MPI_Gather(spec_profiles_pp, spec_ints_pp, MPI_INT, spec_profiles, spec_ints_pp, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(traj_derivs_pp, traj_ints_pp, MPI_DOUBLE, traj_derivs, traj_ints_pp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    
    // Process simulation data on processor 0
    if(id==0){
        
    // Reshape the data into multidimensional arrays
    int spec_profiles_mda[N_traj][N_record+1][n_specs];
    double traj_derivs_mda[N_traj][N_record+1][n_params];
    
    ind_rec_spec = 0;
    ind_rec_rxns = 0;
    
    for(int i=0; i < N_traj; i++){
        for(int j=0; j<N_record+1; j++){
            
            // Reshape species profile data
            for(int k = 0; k < n_specs; k++){
                spec_profiles_mda[i][j][k] = spec_profiles[ind_rec_spec];
                ind_rec_spec += 1;
            }
            
            // Reshape trajectory derivative data
            for(int k = 0; k < n_params; k++){
                traj_derivs_mda[i][j][k] = traj_derivs[ind_rec_rxns];
                ind_rec_rxns += 1;
            }
        }
    }
    
    // Reshape the data into multidimensional arrays
    double spec_profiles_averages[N_record+1][n_specs];
    double traj_deriv_avgs[N_record+1][n_params];
    double sensitivities[N_record+1][n_specs][n_params];
    
    
        
    }
    
    MPI_Finalize ( );     //  Terminate MPI.
    
    if(id==0){
        cout << "Simulation completed successfully" << endl;
    }


  
    /*
    ============ Perform statistical analysis ==============
    */

    // Resize vectors to be able to hold the necessary data
    
    spec_profiles_averages.resize(in_data.N_record);
    for (int i = 0; i < in_data.N_record; i++){
        spec_profiles_averages[i].resize(in_data.n_specs);}

    sensitivities.resize(in_data.n_specs);
    for (int i = 0; i < in_data.n_specs; i++){
        
        sensitivities[i].resize(in_data.N_record);
        for(int j = 0; j < in_data.N_record; j++){
            
            sensitivities[i][j].resize(in_data.n_params);
        }
    }
    
    // Fill in with fake data
    
    for (int i = 0; i < in_data.N_record; ++i){
        for(int j = 0; j < in_data.n_specs; j++){
            spec_profiles_averages[i][j] = 0;
        }
    }

    for (int i = 0; i < in_data.n_specs; ++i){
        for(int j = 0; j < in_data.N_record; j ++){
            for(int k = 0; k < in_data.n_params; k++){
                sensitivities[i][j][k] = 0;
            }
        }
    }
    
    
    // Average species numbers accross trajectories
    for(int i=0; i < N_record+1; i++){
        for(int j=0; j<n_specs; j++){
            
            spec_profiles_averages[i][j] = 0;
            
            for(int rep_num = 0; rep_num < N_traj; rep_num++){
                spec_profiles_averages[i][j] += spec_profiles_mda[rep_num][i][j];
            }
            
            spec_profiles_averages[i][j] = spec_profiles_averages[i][j] / N_traj;
        }
    }
    
    // Average trajectory derivatives accross trajectories
    for(int i=0; i < N_record+1; i++){
        for(int j=0; j < n_params; j++){
            
            traj_deriv_avgs[i][j] = 0;
            
            for(int rep_num = 0; rep_num < N_traj; rep_num++){
                traj_deriv_avgs[i][j] += traj_derivs_mda[rep_num][i][j];
            }
            
            traj_deriv_avgs[i][j] = traj_deriv_avgs[i][j] / N_traj;
        }
    }
    
    // Compute sensitivities as the covariance of species populations and trajectory derivatives
    for(int i=0; i<N_record+1; i++){
        for(int j=0; j<n_specs; j++){
            for(int k = 0; k < n_params; k++){
                
                sensitivities[i][j][k] = 0;
                
                for(int traj_ind = 0; traj_ind < N_traj; traj_ind++){
                    sensitivities[i][j][k] += ( spec_profiles_mda[traj_ind][i][j] - spec_profiles_averages[i][j] ) * ( traj_derivs_mda[traj_ind][i][k]- traj_deriv_avgs[i][k] );
                }
                
                sensitivities[i][j][k] = sensitivities[i][j][k] / (N_traj - 1);
                
            }
            
        }
    }

    // Write output files
    if(id==0){
        write_spec_avg_output();        
        write_sensitivity_output();
    }
    
}