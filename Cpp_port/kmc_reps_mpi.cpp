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


int main ( int argc, char *argv[] );

int main ( int argc, char *argv[] ){
	
/*
============================ Read input file ============================
*/

// Variables required for KMC simulation
int n_rxns;
int n_specs;
int n_params;
double t_final;                         // termination time (s)
int N_record = 100;                     // number of time points to record
int N_traj = 1000;                       // number of replicate trajectories
double eps;                	            // stiffness level, eps << 1
bool write_traj_files = false;           // flag to write data files for each individual trajectory, will take up a lot of space
string* spec_names = NULL;              // names of chemical species
int* N_0 = NULL;                        // initial state
//int** stoich_mat = NULL;                // stoichiometry matrix                      
double* rate_const = NULL;              // rate constants of elementary reactions
string* param_names = NULL;             // names of the parameters

// Two time scale variables
bool two_time_scale = false;            // simulation uses one time scale by default
int* fast_rxns = NULL;                  // indices of the fast reactions
int num_batches = 50;        	        // number of batches to use for microscale steady-state averaging in KMC_TTS
double delta = 0.05;              	    // accuracy level of microscale averaging, delta << 1


int stoich_mat[3][3] = { {-1, 1, 0},
    {1, -1, 0},
    {0, -1, 1}};                       // stoichiometry matrix, only thing we cannot read properly right now


ifstream myfile ("input.txt");
string line;

if (myfile.is_open())
{
    while ( getline (myfile,line) )
    {
        if(line.find_first_not_of(' ') == std::string::npos or line[0] == '#')    // skip empty lines and comments
        {
            continue;
        }
        
        
        if(line == "Number of species"){
            
            getline (myfile,line);
            n_specs = atoi(line.c_str());
            
            spec_names = new string[n_specs];       // allocate memory for the species names
            N_0 = new int[n_specs];                 // allocate memory for the initial state
        
        }else if(line == "Number of reactions"){
            
            getline (myfile,line);
            n_rxns = atoi(line.c_str());
            n_params = n_rxns;
            
            param_names = new string[n_params];       // allocate memory for the species names
            rate_const = new double[n_rxns];
        
        }else if(line == "Species names"){
            
            for(int spec_ind = 0; spec_ind < n_specs; spec_ind++){
                getline (myfile,line);
                spec_names[spec_ind] = line;
            }
            
        }else if(line == "Initial state"){
            
            for(int spec_ind = 0; spec_ind < n_specs; spec_ind++){
                getline (myfile,line);
                N_0[spec_ind] = atoi(line.c_str());
            }
             
        }else if(line == "Final time"){
        
            getline (myfile,line);
            t_final = atof(line.c_str());
            
        //}else if(line == "Reactions"){          // read in stoichiometric matrix
        //    
        //    stoich_mat = new int*[n_rxns];          // allocate stoichiometric matrix
        //    for(int i = 0; i < n_rxns; ++i){
        //        stoich_mat[i] = new int[n_specs];
        //    }
        //        
        //    for(int rxn_ind = 0; rxn_ind < n_rxns; rxn_ind++){
        //        
        //        getline (myfile,line);
        //        string delimiter = " ";
        //        //string token = line.substr(0, line.find(" ")); // token is "scott"
        //        size_t pos = 0;
        //        string token;
        //        
        //        for(int spec_ind = 0; spec_ind < n_specs; spec_ind++){
        //            pos = line.find(delimiter);
        //            token = line.substr(0, pos);
        //            stoich_mat[rxn_ind][spec_ind] = atoi(token.c_str());
        //            line.erase(0, pos + delimiter.length());
        //        } 
        //        
        //    }
            
        }else if(line == "Rate constants"){
            
            for(int rxn_ind = 0; rxn_ind < n_rxns; rxn_ind++){
                getline (myfile,line);
                rate_const[rxn_ind] = atof(line.c_str());
            }
            
        }else if(line == "write trajectory data"){      // Flag to print out data for individual trajectories
            
            write_traj_files = true;
            
        }else if(line == "two time scale"){      // Flag to turn two time scale mode on
            
            two_time_scale = true;
            
        }else if(line == "Parameter names"){
            
            for(int param_ind = 0; param_ind < n_params; param_ind++){
                getline (myfile,line);
                param_names[param_ind] = line;
            }
            
        }else if(line == "Fast reactions"){
            
            // read the fast reactions
            
        }//else{                  // Unrecognized command. Throw an error. Should also tell them the line number
        //    
        //    cout << "Unrecognized command: " << line << endl;
        //    return -1;
        //    
        //}
    }
}else{
    cout << "Unable to open file" << endl;
    return -1;
}

myfile.close();


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

int N[n_specs];                         // species populations
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
    
    // Initialize species populations
    for(int k =0; k < n_specs; k++){
        N[k] = N_0[k];
    }
    
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
                    props[i] = props[i] * pow (N[j], -stoich_mat[i][j]);
                }
            }
            
            // Fill in derivatives for each parameter
            for(int k = 0; k < n_params; k++){
                if(i==k){
                    
                    prop_ders[i][k] = 1.0;
                    
                    for(int j = 0; j < n_specs; j++){
                        if(stoich_mat[i][j] < 0){//if this species is a reactant, use it to compute the rate
                            prop_ders[i][k] = prop_ders[i][k] * pow (N[j], -stoich_mat[i][j]);
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
                    writer_spec << N[j] << "\t";
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
                 spec_profiles_pp[ind_rec_spec] = N[j];
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
            N[j] += stoich_mat[rxn_to_fire_ind][j];
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
                writer_spec << N[j] << "\t";
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
             spec_profiles_pp[ind_rec_spec] = N[j];
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

/*
============ Perform statistical analysis ==============
*/

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

for(int time_ind = 0; time_ind <= N_record; time_ind ++){
    
    writer_spec_avgs << t_rec[time_ind] << "\t";
    
    for(int spec_ind = 0; spec_ind < n_specs; spec_ind++){
        writer_spec_avgs << spec_profiles_averages[time_ind][spec_ind] << "\t"; // print the mean population of this species at this time
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