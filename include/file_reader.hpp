#pragma once

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

    // Methods
    file_reader();                          // empty constructor
    file_reader(string flname);             // constructor which reads the input file

		// Variables required for KMC simulation
		int rand_seed = 1;
    int n_rxns = 0;
		int n_specs = 0;
		int n_params = 0;
		double t_final = 0;                         // termination time (s)
		int N_record = 101;                           // number of time points to record
		int N_traj = 1000;                             // number of replicate trajectories
		double eps = 0;                	            // stiffness level, eps << 1
		bool write_traj_files = false;                  // flag to write data files for each individual trajectory, will take up a lot of space
		vector<string> spec_names;              // names of chemical species
		vector<int> N_0;                        // initial state
		vector< vector<int> > stoich_mat;       // stoichiometric matrix
		vector<double> rate_const;              // rate constants of elementary reactions
		vector<string> param_names;             // names of the parameters
		vector<double> t_rec;                   // vector of sampling times

		// Two time scale variables
		bool two_time_scale = false;                    // simulation uses one time scale by default

    int n_fast_pairs = 0;                       // number of pairs of fast reactions
    vector < vector<int> > fast_pairs;      // indices of fast reaction pairs

		vector<int> fast_rxns;                  // indices of the fast reactions
    vector<int> slow_rxns;                  // indices of the slow reactions
    int n_fast_rxns = 0;                        // number of fast reactions
    int n_slow_rxns = 0;                        // number of slow reactions

    int n_micro_steps = 1000;                      // number of Metropolis steps to use for averaging on the micro scale

		//int num_batches;        	            // number of batches to use for microscale steady-state averaging in KMC_TTS
		//double delta;              	            // accuracy level of microscale averaging, delta << 1





};
