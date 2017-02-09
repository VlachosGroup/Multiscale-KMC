

#ifndef MY_HEADER
#define MY_HEADER

# include <vector>
# include <string>
using namespace std;

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
		vector< vector<int> > stoich_mat;                    
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


class STS_traj {

    private:
    
        static string species_out_flname;
        static string traj_deriv_out_flname;
    
		vector<int> N;                         	// species populations
		vector<double> props;
		vector<double> prop_cum;
		vector<double> asum;
		
		// Used for sensitivity analysis
		vector<double> W;                       // trajectory derivatives
		vector< vector<double> > prop_ders;     // derivative of each propensity with respect to each parameter
		vector<double> prop_ders_sum;           // derivative of the sum of all propensities with respect to each parameter
		double t_trunc;                         // time since the previous KMC step
		
		double t;                               // KMC clock
		double t_prev;                          // KMC time of previous step
		double dt;                              // time step
		double r_rxn_choose;                    // random number between 0 and 1 used to choose which reaction to fire
		double r_timestep;                      // random mumber between 0 and 1 used to choose the time step
		int rxn_to_fire_ind;                    // index of the reaction chosen to fire
		int ind_rec;                            // time point
		int ind_rec_spec;                   	// index to record species population data for all trajectories running on the current processor
		int ind_rec_rxns;                   	// index to record trajectory derivative data for all trajectories running on the current processor
	
		
	
    public:
    
        file_reader in_data;            		// input data from the input file
	
        // Data to record
		// species vs. time
		// trajectory derivative vs. time
    
        STS_traj();                          	// empty constructor
		void simulate(int);								// perform KMC simulation
		void record_data();							// record data at the current KMC time
		
	
};


class TTS_traj : public STS_traj {

    private:
    
		// need an extra 3-d vector of microscale sensitivities
	
    public:

		TTS_traj();                          // empty constructor
};

#endif