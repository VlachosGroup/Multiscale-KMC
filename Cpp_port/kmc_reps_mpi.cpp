# include <iostream>    // using IO functions
# include <string>      // using string
# include <vector>
# include <fstream>
# include <random>
# include <math.h>
# include <mpi.h>
using namespace std;


int main ( int argc, char *argv[] );

int main ( int argc, char *argv[] ){
	
/*
============================ User input parameters ============================
*/

// Need to read this data from an input file

// Specify the reaction system - used for all models
string spec_names[3] = {"A","B","C"};                   // names of chemical species
int N_0[3] = {100, 0, 0};                             // initial state
int stoich_mat[3][3] = { {-1, 1, 0},
    {1, -1, 0},
    {0, -1, 1}};                       // stoichiometry matrix                      
double rate_const[3] = {20.0, 30.0, 2.0};                            // rate constants of elementary reactions
string param_names[3] = {"k1","k2","k3"};                  // names of the parameters
double t_final = 10;                         // termination time (s)
int N_record = 100;                           // number of time points to record
int fast_rxns[2] = {1,2};                       // indices of the fast reactions
int N_traj = 1000;                               // number of replicate trajectories

// Only for STS ODE and KMC
double eps = 0.01;                	// stiffness level, eps << 1

// Only for TTS KMC
int num_batches = 20;        	// number of batches to use for microscale steady-state averaging in KMC_TTS
double delta = 0.1;              	// accuracy level of microscale averaging, delta << 1

// Should be able to automate getting this from the input
int n_rxns = 3;
int n_specs = 3;
int n_params = n_rxns;

bool write_traj_files = false;      // set to true if the user wants data files for each individual trajectory, will take up a lot of space

/*
============================ Set up parallelization variables ============================
*/

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
============================ Declare variables ============================
*/

double props[n_rxns];
double prop_cum[n_rxns];
double asum;

// Used for sensitivity analysis
double W[n_params];                     // trajectory derivatives
double prop_ders[n_rxns][n_params];     // derivative of each propensity with respect to each parameter
double prop_ders_sum[n_params];         // derivative of the sum of all propensities with respect to each parameter
double t_trunc;                         // time since the previous KMC step

double t;                               // KMC clock
double t_prev;                          // KMC time of previous step
double dt;                              // time step
double r_rxn_choose;                    // random number between 0 and 1 used to choose which reaction to fire
double r_timestep;                      // random mumber between 0 and 1 used to choose the time step
int rxn_to_fire_ind;                    // index of the reaction chosen to fire
double t_rec[N_record+1];               // array of times at which to record data
int ind_rec;                            // time point
int ind_rec_spec = 0;                   // index to record species population data for all trajectories running on the current processor
int ind_rec_rxns = 0;                   // index to record trajectory derivative data for all trajectories running on the current processor

// Fill in recording times
for(int i = 0; i <= N_record; i++){
    t_rec[i] = (double) i / N_record * t_final;
}

for(int rep=0; rep < Npp; rep++){
    
    srand(sub_seeds[rep]);      // Set the random seed
    
    /*
    ============================ Initialize variables ============================
    */
    
    t = 0;
    t_prev = 0;
    ind_rec = 0;
    
    // Initialize trajectory derivatives
    for(int k =0; k < n_params; k++){
        W[k] = 0;
    }
    
    // Open files to record data for the trajectory
    ofstream writer_spec;
    ofstream writer_SA;
    
    if(write_traj_files){
        
        // Create a folder for all trajectory data files
        // Name the file specific to this trajectory. Put random seed in the file name.
        
        writer_spec.open("kmc_specs.txt");
        if(! writer_spec){  
            cout << "Error opening file" << endl;
            return -1;
        }
        
        // Write header    
        writer_spec << "Time \t";
        for(int j = 0; j < n_specs; j++){
            writer_spec << spec_names[j] << "\t";
        }
        writer_spec << endl;
        
        writer_SA.open("kmc_SA.txt");
        if(! writer_SA){  
            cout << "Error opening file" << endl;
            return -1;
        }
        
        // Write header
        writer_SA << "Time \t";
        for(int i = 0; i < n_params; i++){
            writer_SA << param_names[i] << "\t";
        }
        writer_SA << endl; 
    }        
    
    // Start KMC loop
    while(t < t_final){
        
        /*
        ============================ Compute quantities for current step ============================
        */
        
        // Generate random numbers
        r_rxn_choose = ((double) rand() / (RAND_MAX));
        r_timestep = ((double) rand() / (RAND_MAX));
        
        // Compute reaction propensities
        for(int i = 0; i < n_rxns; i++){
            
            // Use mass action kinetics equation, a = k * [A]^ma * [B]^mb * ...
            props[i] = rate_const[i];
            for(int j = 0; j < n_specs; j++){
                if(stoich_mat[i][j] < 0){//if this species is a reactant, use it to compute the rate
                    props[i] = props[i] * pow (N_0[j], -stoich_mat[i][j]);
                }
            }
            
            // Fill in derivatives for each parameter
            for(int k = 0; k < n_params; k++){
                if(i==k){
                    
                    prop_ders[i][k] = 1.0;
                    
                    for(int j = 0; j < n_specs; j++){
                        if(stoich_mat[i][j] < 0){//if this species is a reactant, use it to compute the rate
                            prop_ders[i][k] = prop_ders[i][k] * pow (N_0[j], -stoich_mat[i][j]);
                        }
                    }
                    
                }else{
                    prop_ders[i][k] = 0;
                }
            }
        }
        
        // Add derivatives for each parameter
        for(int k=0; k < n_params; k++){
            prop_ders_sum[k] = 0;
            for(int i=0; i < n_rxns; i++){
                prop_ders_sum[k] += prop_ders[i][k];
            }
        }
        
        // Sum the propensities
        asum = 0;
        for(int i = 0; i < n_rxns; i++){
            asum += props[i];
        }
        
        // If all propensities are 0, then exit the while loop
        if(asum == 0){
            break;
        }
        
        // Choose which reaction to fire
        for(int i = 0; i < n_rxns; i++){
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
        ============================ Record data ============================
        */
        
        // Record the current state as long as time >= t_sample
        while(t >= t_rec[ind_rec]){
            
            
            t_trunc = t_rec[ind_rec] - t_prev;
            
            // Record data in output files
            if(write_traj_files){
                
                // Record the species numbers
                writer_spec << t_rec[ind_rec] << "\t";
                for(int j = 0; j < n_specs; j++){
                    writer_spec << N_0[j] << "\t";
                }
                writer_spec << endl;
            
                // Record trajectory derivatives
                writer_SA << t_rec[ind_rec] << "\t";
                for(int k = 0; k < n_params; k++){
                    writer_SA << W[k] - prop_ders_sum[k] * t_trunc << "\t";
                }
                writer_SA << endl;
            }
            
            // Record data in variables
            
            // Record species populations
            for(int j = 0; j < n_specs; j++){
                 spec_profiles_pp[ind_rec_spec] = N_0[j];
                 ind_rec_spec += 1;
                }
            
            // Record trajectory derivatives
            for(int k = 0; k < n_params; k++){
                traj_derivs_pp[ind_rec_rxns] = W[k] - prop_ders_sum[k] * t_trunc;
                ind_rec_rxns += 1;
                }
            
            
            ind_rec += 1;
            
        }
        
        /*
        ============ Fire reaction and update system ==============
        */
        
        // Update species populations
        for(int j = 0; j < n_specs; j++){
            N_0[j] += stoich_mat[rxn_to_fire_ind][j];
        }
        
        // Update clock
        t_prev = t;
        t = t + dt;
        
        // Update trajectory derivatives
        for(int k=0; k < n_params; k++){
            
            W[k] += prop_ders[rxn_to_fire_ind][k] / props[rxn_to_fire_ind];         // contribution from reaction that fires            
            W[k] -= prop_ders_sum[k] * dt;       // contribution from time step
        }

    }
    
    // Fill in recording times that were missed
    while(ind_rec <= N_record){
        
        t_trunc = t_rec[ind_rec] - t_prev;
            
        // Record data in output files
        if(write_traj_files){
            
            // Record the species numbers
            writer_spec << t_rec[ind_rec] << "\t";
            for(int j = 0; j < n_specs; j++){
                writer_spec << N_0[j] << "\t";
            }
            writer_spec << endl;
        
            // Record trajectory derivatives
            writer_SA << t_rec[ind_rec] << "\t";
            for(int k = 0; k < n_params; k++){
                writer_SA << W[k] - prop_ders_sum[k] * t_trunc << "\t";
            }
            writer_SA << endl;
        }
        
        // Record data in variables
        
        // Record species populations
        for(int j = 0; j < n_specs; j++){
             spec_profiles_pp[ind_rec_spec] = N_0[j];
             ind_rec_spec += 1;
            }
        
        // Record trajectory derivatives
        for(int k = 0; k < n_params; k++){
            traj_derivs_pp[ind_rec_rxns] = W[k] - prop_ders_sum[k] * t_trunc;
            ind_rec_rxns += 1;
            }
        
        
        ind_rec += 1;
            
        }
    
    if(write_traj_files){
        writer_spec.close();
        writer_SA.close();
    }
    

}

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
    
    for(int i=0; i < N_traj; i++){
        for(int j=0; j<N_record+1; j++){
            
            // Reshape species profile data
            for(int k = 0; k < n_specs; k++){
                spec_profiles_mda[i][j][k] = 0;         // PLACEHOLDER VALUE
            }
            
            // Reshape trajectory derivative data
            for(int k = 0; k < n_params; k++){
                traj_derivs_mda[i][j][k] = 0;           // PLACEHOLDER VALUE
            }
        }
    }
    
    // Reshape the data into multidimensional arrays
    double spec_profiles_averages[N_record+1][n_specs];
    double sensitivities[N_record+1][n_specs][n_params];
    
    /*
    ============ Perform statistical analysis ==============
    */

    for(int i=0; i<N_record+1; i++){
        for(int j=0; j<n_specs; j++){
            spec_profiles_averages[i][j] = 0;           // PLACEHOLDER VALUE
        }
    }
    
    for(int i=0; i<N_record+1; i++){
        for(int j=0; j<n_specs; j++){
            for(int k = 0; k < n_params; k++){
                sensitivities[i][j][k] = 0;           // PLACEHOLDER VALUE
            }
            
        }
    }
    
    

    /*
    ============ Print species population averages into an output file ==============
    */
    
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
    
    for(int time_ind = 0; time_ind < N_record+1; time_ind ++){
        
        writer_spec_avgs << t_rec[time_ind] << "\t";
        
        for(int spec_ind = 0; spec_ind < n_specs; spec_ind++){
            writer_spec_avgs << spec_profiles_averages[time_ind][spec_ind] << "\t";
        }
        
        writer_spec_avgs << endl;
    }
    
    writer_spec_avgs.close();
    
    /*
    ============ Print sensitivities into an output file ==============
    */
    
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

MPI_Finalize ( );     //  Terminate MPI.

if(id==0){
    cout << "Simulation completed successfully" << endl;
}

return 0;
}