# include <iostream>    // using IO functions
# include <string>      // using string
# include <vector>
# include <fstream>
# include <random>
# include <math.h>
//# include <mpi.h>

# include <typeinfo>
# include <sstream>
# include <stdio.h>      /* printf, fgets */
# include <stdlib.h>     /* atoi */
# include "myHeader.h"
using namespace std;


class Traj_stats_STS {
	
	private:
	
    static string species_avgs_out_flname;                      // "spec_avgs.out"      
    static string SA_out_flname;                                // "SA.out"
    
    vector< vector< vector<double> > > spec_profiles;           // species data from each trajectory
    vector< vector< vector<double> > > traj_derivs;             // trajectory derivative data from each trajectory
    
    vector< vector<double> > spec_profiles_averages;            // species averages
    vector< vector< vector<double> > > sensitivities;           // sensities
    
	public:
    
    file_reader in_data;            // input data from the input file
    
    // Empty class constructor
    Traj_stats_STS() : in_data(){}
    
    /*
    ============ Resize vectors to be able to hold the necessary data =============
    */
    void allocate_data(){
        
        spec_profiles.resize(in_data.N_traj);
        for (int i = 0; i < in_data.N_traj; i++){
            
            spec_profiles[i].resize(in_data.N_record);
            for(int j = 0; j < in_data.N_record; j ++){
                
                spec_profiles[i][j].resize(in_data.n_specs);
            }
        }
        
        traj_derivs.resize(in_data.N_traj);
        for (int i = 0; i < in_data.N_traj; i++){
            
            traj_derivs[i].resize(in_data.N_record);
            for(int j = 0; j < in_data.N_record; j ++){
                
                traj_derivs[i][j].resize(in_data.n_params);
            }
        }
        
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
        
    }
    
    /*
    ============================ Run simulations to gather data ============================
    */
    
    void run_simulations(){
        
        // Fill in fake data
        for (int i = 0; i < in_data.N_traj; ++i){
            for(int j = 0; j < in_data.N_record; j ++){
                for(int k = 0; k < in_data.n_specs; k++){
                    spec_profiles[i][j][k] = 0;
                }
            }
        }
        
        for (int i = 0; i < in_data.N_traj; ++i){
            for(int j = 0; j < in_data.N_record; j ++){
                for(int k = 0; k < in_data.n_params; k++){
                    traj_derivs[i][j][k] = 0;
                }
            }
        }
        
        STS_traj run;
        run.in_data = in_data;      // Copy input file data to the trajectory object
        run.simulate(12345);
        
//        // Set up MPI variables
//        
//        int id;
//        int ierr;
//        int p;
//        int Npp;
//        
//        ierr = MPI_Init ( &argc, &argv );                 //  Initialize MPI.
//        ierr = MPI_Comm_size ( MPI_COMM_WORLD, &p );      //  Get the number of processes.  
//        ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id );     //  Get the individual process ID.
//        
//        
//        
//        Npp = N_traj / p + 1;
//        N_traj = Npp * p;       // round up to the nearest multiple of the number of processors
//        
//        if(id==0){
//            cout << N_traj << " replicate trajectories will be used." << endl;
//        }
//        
//        
//        int seed_list[N_traj];
//        for(int i = 0; i < N_traj; i++){
//            seed_list[i] = i+1;
//        }
//        
//        
//        // For each processor
//        int spec_ints_pp = Npp * (N_record+1) * n_specs;
//        int traj_ints_pp = Npp * (N_record+1) * n_params;
//        int spec_profiles_pp[spec_ints_pp];
//        double traj_derivs_pp[traj_ints_pp];
//        // Fill these up with data as the simulation progresses
//        
//        // Scatter the seed list
//        int sub_seeds[Npp];
//        MPI_Scatter(seed_list, Npp, MPI_INT, sub_seeds,
//                Npp, MPI_INT, 0, MPI_COMM_WORLD);
//        
//
//        /*
//        ============ run KMC simulations to gather data ==============
//        */       
//        
//        
//        /*
//        ============ Collect data from the replicate simulations ==============
//        */
//        
//        // For final analysis, gather data into the big arrays
//        int spec_profiles[N_traj * (N_record+1) * n_specs];
//        double traj_derivs[N_traj * (N_record+1) * n_params];
//        
//        MPI_Gather(spec_profiles_pp, spec_ints_pp, MPI_INT, spec_profiles, spec_ints_pp, MPI_INT, 0, MPI_COMM_WORLD);
//        MPI_Gather(traj_derivs_pp, traj_ints_pp, MPI_DOUBLE, traj_derivs, traj_ints_pp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//        
//        
//        // Process simulation data on processor 0
//        if(id==0){
//            
//        // Reshape the data into multidimensional arrays
//        int spec_profiles_mda[N_traj][N_record+1][n_specs];
//        double traj_derivs_mda[N_traj][N_record+1][n_params];
//        
//        ind_rec_spec = 0;
//        ind_rec_rxns = 0;
//        
//        for(int i=0; i < N_traj; i++){
//            for(int j=0; j<N_record+1; j++){
//                
//                // Reshape species profile data
//                for(int k = 0; k < n_specs; k++){
//                    spec_profiles_mda[i][j][k] = spec_profiles[ind_rec_spec];
//                    ind_rec_spec += 1;
//                }
//                
//                // Reshape trajectory derivative data
//                for(int k = 0; k < n_params; k++){
//                    traj_derivs_mda[i][j][k] = traj_derivs[ind_rec_rxns];
//                    ind_rec_rxns += 1;
//                }
//            }
//        }
//        
//        // Reshape the data into multidimensional arrays
//        double spec_profiles_averages[N_record+1][n_specs];
//        double traj_deriv_avgs[N_record+1][n_params];
//        double sensitivities[N_record+1][n_specs][n_params];
//        
//        
//            
//        }
//        
//        MPI_Finalize ( );     //  Terminate MPI.
//        
//        if(id==0){
//            cout << "Simulation completed successfully" << endl;
//        }

    }
    
    /*
    ============ Perform statistical analysis ==============
    */
    void process_statistics(){
        
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
        
        
//        // Average species numbers accross trajectories
//        for(int i=0; i < N_record+1; i++){
//            for(int j=0; j<n_specs; j++){
//                
//                spec_profiles_averages[i][j] = 0;
//                
//                for(int rep_num = 0; rep_num < N_traj; rep_num++){
//                    spec_profiles_averages[i][j] += spec_profiles_mda[rep_num][i][j];
//                }
//                
//                spec_profiles_averages[i][j] = spec_profiles_averages[i][j] / N_traj;
//            }
//        }
//        
//        // Average trajectory derivatives accross trajectories
//        for(int i=0; i < N_record+1; i++){
//            for(int j=0; j < n_params; j++){
//                
//                traj_deriv_avgs[i][j] = 0;
//                
//                for(int rep_num = 0; rep_num < N_traj; rep_num++){
//                    traj_deriv_avgs[i][j] += traj_derivs_mda[rep_num][i][j];
//                }
//                
//                traj_deriv_avgs[i][j] = traj_deriv_avgs[i][j] / N_traj;
//            }
//        }
//        
//        // Compute sensitivities as the covariance of species populations and trajectory derivatives
//        for(int i=0; i<N_record+1; i++){
//            for(int j=0; j<n_specs; j++){
//                for(int k = 0; k < n_params; k++){
//                    
//                    sensitivities[i][j][k] = 0;
//                    
//                    for(int traj_ind = 0; traj_ind < N_traj; traj_ind++){
//                        sensitivities[i][j][k] += ( spec_profiles_mda[traj_ind][i][j] - spec_profiles_averages[i][j] ) * ( traj_derivs_mda[traj_ind][i][k]- traj_deriv_avgs[i][k] );
//                    }
//                    
//                    sensitivities[i][j][k] = sensitivities[i][j][k] / (N_traj - 1);
//                    
//                }
//                
//            }
//        }

    }
    
    
    /* 
    =========== Print species population averages into an output file ===========
    */
    void write_spec_avg_output(){

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
    void write_sensitivity_output(){
        
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
                    writer_sensitivities << sensitivities[i][time_ind][param_ind] << "\t";
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
    
};

string Traj_stats_STS :: species_avgs_out_flname = "spec_avgs.out";
string Traj_stats_STS :: SA_out_flname = "SA.out";


class Traj_stats_TTS : public Traj_stats_STS{
    
    private:
    
    double microscale_sensitivity_contributions[2][2][2];
    
    public:
    
    Traj_stats_TTS() : Traj_stats_STS(){
        
    }
    
    void add_microscale_sensitivities(){
        // add microscale contributions to the sensitivities
    }
};


// Main function
int main() {
	
	file_reader fr("network.in");               // Read input file
        
        
    if(! fr.two_time_scale){                    // single time scale
            
        Traj_stats_STS STS_sim; 
        STS_sim.in_data = fr;                   // assign input
        STS_sim.allocate_data();
            
        STS_sim.run_simulations();              // run KMC simulations
        STS_sim.process_statistics();
        
        STS_sim.write_spec_avg_output();        // Write output files
        STS_sim.write_sensitivity_output();
        
    }else{                                      // two time scale version
        
        Traj_stats_TTS TTS_sim;
        TTS_sim.in_data = fr;
        TTS_sim.allocate_data();
        
        TTS_sim.run_simulations();
        TTS_sim.process_statistics();
        TTS_sim.add_microscale_sensitivities();
        
        TTS_sim.write_spec_avg_output();
        TTS_sim.write_sensitivity_output();
        
    }
    
    cout << "Simulation complete." << endl;

	return 0;
}