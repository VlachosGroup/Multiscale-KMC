#ifndef MY_HEADER
#define MY_HEADER

# include <vector>
# include <string>
# include <iostream>    // using IO functions
# include <string>      // using string
# include <fstream>
# include <sstream>
using namespace std;


/*
============================ Class to read input file ============================
*/

class file_reader {

    private:
    
    public:

		// Variables required for KMC simulation
		int n_rxns;
		int n_specs;
		int n_params;
		double t_final;                         // termination time (s)
		int N_record;                           // number of time points to record
		int N_traj;                             // number of replicate trajectories
		double eps;                	            // stiffness level, eps << 1
		bool write_traj_files;                  // flag to write data files for each individual trajectory, will take up a lot of space
		vector<string> spec_names;              // names of chemical species
		vector<int> N_0;                        // initial state
		vector< vector<int> > stoich_mat;       // stoichiometric matrix
		vector<double> rate_const;              // rate constants of elementary reactions
		vector<string> param_names;             // names of the parameters
		vector<double> t_rec;                   // vector of sampling times
        
		// Two time scale variables
		bool two_time_scale;                    // simulation uses one time scale by default
        
        int n_fast_pairs;                       // number of pairs of fast reactions
        vector < vector<int> > fast_pairs;      // indices of fast reaction pairs
        
		vector<int> fast_rxns;                  // indices of the fast reactions
        vector<int> slow_rxns;                  // indices of the slow reactions
        int n_fast_rxns;                        // number of fast reactions
        int n_slow_rxns;                        // number of slow reactions
        
        int n_micro_steps;                      // number of Metropolis steps to use for averaging on the micro scale
        
		//int num_batches;        	            // number of batches to use for microscale steady-state averaging in KMC_TTS
		//double delta;              	            // accuracy level of microscale averaging, delta << 1
	
    
        // Methods
		file_reader(string flname);             // constructor which reads the input file
		file_reader();                          // empty constructor
        
};


/*
============================ Class to handle many trajectories ============================
*/

class Traj_stats {
	
	private:
	
        static string species_avgs_out_flname;                      // "spec_avgs.out"
        static string SA_out_flname;                                // "SA.out"
        
        vector< vector<double> > spec_profiles_averages;            // species averages
        vector< vector<double> > traj_deriv_avgs;                   // trajectory derivative averages
        vector< vector< vector<double> > > sensitivities;           // sensities
        vector< vector< vector<double> > > microscale_sensitivity_contributions;
    
	public:
    
        file_reader in_data;            // input data from the input file
        
        Traj_stats();               // Empty class constructor
        void initialize_stats();
        void run_simulations();
        void finalize_stats();
        void write_spec_avg_output();
        void write_sensitivity_output();
    
};



/*
============================ Class to handle one STS trajectory ============================
*/


class KMC_traj {

    private:
    
    public:
    
        // File writer objects for output files
        static string species_out_flname;
        static string traj_deriv_out_flname;
        
        // Basic KMC variables
        vector <int> N;
        double t_trunc;                                         // time since the previous KMC step
        double t;                                               // KMC clock
        double t_prev;                                          // KMC time of previous step
        double dt;                                              // time step
        int rxn_to_fire_ind;                                    // index of the reaction chosen to fire
        int ind_rec;                                            // time point
        vector <double> props;
        
        // Used for sensitivity analysis
        vector <double> W;                                      // trajectory derivatives
        vector < vector <double> > prop_ders;     // derivative of each propensity with respect to each parameter
        vector <double> prop_ders_sum;
        
    
        // For file writing
        ofstream writer_spec;
        ofstream writer_SA;
        
        
        void record_stats();
    
    
    
        file_reader in_data;            		// input data from the input file
	
        // Data to record
		vector< vector<double> > spec_profile;                      // species vs. time
		vector< vector<double> > traj_deriv_profile;                // trajectory derivative vs. time
    
        KMC_traj();                          	             // empty constructor
		void simulate(int);								// perform STS KMC simulation
        void initialize_sim(int);
	
};

/*
============================ Class to handle one TTS trajectory ============================
*/

class KMC_traj_TTS : public KMC_traj {

    private:
    
        vector <double> N_micro_avg;
        vector < vector <double> > micro_scale_sens;
        
        vector < vector <double> > prop_ders_direct;        // Direct averaging of derviativeson the microscale
        vector < vector <double> > prop_ders_indirect;      // Indirect averaging on the microscale
        
        void record_stats_TTS();
		void simulate_micro();								// perform microscale KMC simulation for fast reactions
	
    public:

        vector< vector< vector<double> > > micro_scale_sens_profile;         // 3-D vector of microscale sensitivities
    
		KMC_traj_TTS();                          // empty constructor
        void simulate_TTS(int);
        
};

#endif