# include <iostream>    // using IO functions
# include <string>      // using string
# include <vector>
# include <fstream>
# include <random>
# include <math.h>
using namespace std;


int main() {
	
    /*
    ============================ Input parameters ============================
    */
    
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
	
	// Only for STS ODE and KMC
	double eps = 0.01;                	// stiffness level, eps << 1
	
	// Only for TTS KMC
	int num_batches = 20;        	// number of batches to use for microscale steady-state averaging in KMC_TTS
	double delta = 0.1;              	// accuracy level of microscale averaging, delta << 1
    
    /*
    ============================ Calculate a few additional things we need ============================
    */

	int n_rxns = 3;
    int n_specs = 3;
    int n_params = n_rxns;
    
    double props[n_rxns];
    double prop_cum[n_rxns];
    double asum;
    
    // Used for sensitivity analysis
    double W[n_params];
    double prop_ders[n_rxns][n_params];
    double prop_ders_sum[n_params];
    double t_trunc;
    
    double t = 0;
    double t_prev = 0;
    double dt;
    double r_rxn_choose;
    double r_timestep;
    int rxn_to_fire_ind;
    double t_rec[N_record];
    int ind_rec = 0;
    
    
    // Initialize trajectory derivatives
    for(int k =0; k < n_params; k++){
        W[k] = 0;
    }
    
    // Fill in recording times
    for(int i = 0; i <= N_record; i++){
        t_rec[i] = (double) i / N_record * t_final;
    }
    
    ofstream writer_spec("kmc_specs.txt");
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

    
    ofstream writer_SA("kmc_SA.txt");
    if(! writer_SA){  
        cout << "Error opening file" << endl;
        return -1;
    }
    
    // Write header
    writer_SA << "Time \t";
    for(int i = 0; i < n_specs; i++){
        writer_SA << param_names[i] << "\t";
    }
    writer_SA << endl;    
    
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
        ============================Record data ============================
        */
        
        // Record the current state as long as time >= t_sample
        while(t >= t_rec[ind_rec]){
            
            // Record the species numbers
            writer_spec << t_rec[ind_rec] << "\t";
            for(int j = 0; j < n_specs; j++){
                writer_spec << N_0[j] << "\t";
            }
            writer_spec << endl;
            
            // Record trajectory derivatives
            t_trunc = t_rec[ind_rec] - t_prev;
            
            writer_SA << t_rec[ind_rec] << "\t";
            for(int k = 0; k < n_params; k++){
                writer_SA << W[k] - prop_ders_sum[k] * t_trunc << "\t";
            }
            writer_SA << endl;
            
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
            
        // Record the species numbers
        writer_spec << t_rec[ind_rec] << "\t";
        for(int j = 0; j < n_specs; j++){
            writer_spec << N_0[j] << "\t";
        }
        writer_spec << endl;
        
        // Record trajectory derivatives
        t_trunc = t_rec[ind_rec] - t_prev;
        
        writer_SA << t_rec[ind_rec] << "\t";
        for(int k = 0; k < n_params; k++){
            writer_SA << W[k] - prop_ders_sum[k] * t_trunc << "\t";
        }
        writer_SA << endl;
        
        ind_rec += 1;
            
        }
    
    writer_spec.close();
    writer_SA.close();
    
	return 0;
}