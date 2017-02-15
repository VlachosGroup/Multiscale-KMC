#ifndef MY_HEADER
#define MY_HEADER

# include <vector>
# include <string>
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
		vector<int> fast_rxns;                  // indices of the fast reactions
		int num_batches;        	            // number of batches to use for microscale steady-state averaging in KMC_TTS
		double delta;              	            // accuracy level of microscale averaging, delta << 1
	
	
		file_reader(string flname);             // constructor which reads the input file
		file_reader();                          // empty constructor
};


/*
============================ Class to handle many STS trajectories ============================
*/

class Traj_stats_STS {
	
	private:
	
        static string species_avgs_out_flname;                      // "spec_avgs.out"
        static string SA_out_flname;                                // "SA.out"
        
        vector< vector<double> > spec_profiles_averages;            // species averages
        vector< vector<double> > traj_deriv_avgs;                   // trajectory derivative averages
        vector< vector< vector<double> > > sensitivities;           // sensities
    
	public:
    
        file_reader in_data;            // input data from the input file
        
        Traj_stats_STS();               // Empty class constructor
        void initialize_stats();
        void run_simulations();
        void finalize_stats();
        void write_spec_avg_output();
        void write_sensitivity_output();
    
};


/*
============================ Class to handle many TTS trajectories ============================
*/

class Traj_stats_TTS : public Traj_stats_STS{
    
    private:
    
    vector< vector< vector<double> > > microscale_sensitivity_contributions;
    
    public:
    
    Traj_stats_TTS() : Traj_stats_STS(){
        
    }
    
    void add_microscale_sensitivities(){
        // add microscale contributions to the sensitivities
    }
};


/*
============================ Class to handle one STS trajectory ============================
*/


class STS_traj {

    private:
    
        static string species_out_flname;
        static string traj_deriv_out_flname;
    
    public:
    
        file_reader in_data;            		// input data from the input file
	
        // Data to record
		vector< vector<double> > spec_profile;                      // species vs. time
		vector< vector<double> > traj_deriv_profile;                // trajectory derivative vs. time
    
        STS_traj();                          	// empty constructor
		void simulate(int);								// perform KMC simulation
		
	
};


/*
============================ Class to handle one TTS trajectory ============================
*/

class TTS_traj : public STS_traj {

    private:
    
		vector< vector< vector<double> > > micro_scale_sens;         // 3-D vector of microscale sensitivities
	
    public:

		TTS_traj();                          // empty constructor
};

#endif