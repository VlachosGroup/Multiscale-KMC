#pragma once

# include <vector>
# include <string>
# include <iostream>    // using IO functions
# include <string>      // using string
# include <fstream>
# include <sstream>

# include "KMC_traj.hpp"

using namespace std;

/*
============================ Class to handle one TTS trajectory ============================
*/

class KMC_traj_TTS : public KMC_traj {

    private:

        vector <double> N_micro_avg;
        vector < vector <double> > micro_scale_sens;

        void record_stats_TTS();
		void simulate_micro();								// perform microscale KMC simulation for fast reactions

        // Microscale vectors
        vector <double> dEdth;
        vector <double> dEdth_avg;
        vector <double> rev_prop_ders;
        vector <vector<int>> N_rec;
        vector <double> q_cum;
        vector <int> N_candidate;
        vector <double> micro_props;
        vector <vector<double>> micro_prop_ders;
        vector <vector<double>> prop_ders_direct;           // Direct averaging of derviativeson the microscale
        vector <vector<double>> prop_ders_indirect;         // Indirect averaging on the microscale
        vector <double> slow_props;
        vector <double> slow_props_cum;

    public:

        vector< vector< vector<double> > > micro_scale_sens_profile;         // 3-D vector of microscale sensitivities

		KMC_traj_TTS();                          // empty constructor
        void simulate(int);
        double get_micro_scale_sens_profile(int, int, int);
};
