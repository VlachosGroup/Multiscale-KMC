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

		// Variables required for KMC simulation
		int rand_seed;
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
