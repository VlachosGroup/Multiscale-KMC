

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
	
	// Two time scale variables
	bool two_time_scale;                    // simulation uses one time scale by default
	vector<int> fast_rxns;                  // indices of the fast reactions
	int num_batches;        	            // number of batches to use for microscale steady-state averaging in KMC_TTS
	double delta;              	            // accuracy level of microscale averaging, delta << 1


	file_reader(string flname);             // constructor which reads the input file
    file_reader();                          // empty constructor
};

#endif