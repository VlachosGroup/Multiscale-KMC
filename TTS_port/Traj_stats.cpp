# include <iostream>    // using IO functions
# include <string>      // using string
# include <vector>
# include <fstream>
# include <random>
# include <math.h>
# include <mpi.h>

# include <typeinfo>
# include <sstream>
# include <stdio.h>      /* printf, fgets */
# include <stdlib.h>     /* atoi */
using namespace std;


class Traj_stats_STS {
	
	private:
	
    file_reader in_data;            // input data from the input file
    
    int spec_profiles_mda[][][]      // species data from each trajectory
    double traj_derivs_mda[][][]     // trajectory derivative data from each trajectory
    
    double spec_profiles_averages[time_ind][spec_ind]   // species averages
    double sensitivities[i][j][k]                       // sensities
    
	public:
    
    /*
    ============================ Class constructor ============================
    */
    file_reader_STS()(string flname){
        // construct this class
    }
    
    /*
    ============================ Run simulations to gather data ============================
    */
    
    void run_simulations(){
        
        // Set up MPI variables
        
        int id;
        int ierr;
        int p;
        int Npp;
        
        ierr = MPI_Init ( &argc, &argv );                 //  Initialize MPI.
        ierr = MPI_Comm_size ( MPI_COMM_WORLD, &p );      //  Get the number of processes.  
        ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id );     //  Get the individual process ID.
        
        
        
        Npp = N_traj / p;
        N_traj = Npp * p;       // round down to the nearest multiple of the number of processors
        
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
        

        /*
        ============ run KMC simulations to gather data ==============
        */       
        
        
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
    }
    
    /*
    ============ Perform statistical analysis ==============
    */
    void process_statistics(){
        
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
    }
    
    
    /* 
    =========== Print species population averages into an output file ===========
    */
    void write_spec_avg_output(){

        ofstream writer_spec_avgs;
        
        writer_spec_avgs.open("spec_avgs.txt");
        if(! writer_spec_avgs){  
            cout << "Error opening file" << endl;
            return -1;
        }
        
        // Write header
        writer_spec_avgs << "Species mean populations versus time" << endl;
        writer_spec_avgs << endl;
        writer_spec_avgs << "Time \t";
        for(int j = 0; j < n_specs; j++){
            writer_spec_avgs << spec_names[j] << "\t";
        }
        writer_spec_avgs << endl;
        
        for(int time_ind = 0; time_ind <= N_record; time_ind ++){
            
            writer_spec_avgs << t_rec[time_ind] << "\t";
            
            for(int spec_ind = 0; spec_ind < n_specs; spec_ind++){
                writer_spec_avgs << spec_profiles_averages[time_ind][spec_ind] << "\t"; // print the mean population of this species at this time
            }
            
            writer_spec_avgs << endl;
        }
        
        writer_spec_avgs.close();
    }
            
    
    /*
    ===========Print sensitivities into an output file ===========
    */
    void write_spec_avg_output(){
        
        ofstream writer_sensitivities;

        writer_sensitivities.open("sensitivities.txt");
        if(! writer_sensitivities){  
            cout << "Error opening file" << endl;
            return -1;
        }
        
        // Write header
        writer_sensitivities << "Sensitivities for each species and rate constant versus time" << endl << endl;
        
        // Loop over species
        for(int i=0; i < n_specs; i++){
            
            writer_sensitivities << spec_names[i] << endl << endl;
            
            // Write header for this species
            writer_sensitivities << "Time \t";
            for(int k = 0; k < n_params; k++){
                writer_sensitivities << param_names[k] << "\t";
            }
            writer_sensitivities << endl;
            
            // Loop over parameters to print out 
            for(int time_ind = 0; time_ind < N_record+1; time_ind ++){
                
                writer_sensitivities << t_rec[time_ind] << "\t";
                
                for(int param_ind = 0; param_ind < n_params; param_ind++){
                    writer_sensitivities << sensitivities[time_ind][i][param_ind] << "\t";
                }
                
                writer_sensitivities << endl;
            }
            
            // Put a barrier in between species data
            if(i < n_specs - 1){
                writer_sensitivities << endl;
                writer_sensitivities << "=======================================" << endl;
                writer_sensitivities << endl;
            }
            
        }
        
        writer_sensitivities.close();
    }
    
};



class Traj_stats_TTS : public Traj_stats_STS{
    
    private:
    
    double microscale_sensitivity_contributions[][][];
    
    public:
    
    Traj_stats_TTS(){
        
    }
    
    void add_microscale_sensitivities(){
        // add microscale contributions to the sensitivities
    }
};