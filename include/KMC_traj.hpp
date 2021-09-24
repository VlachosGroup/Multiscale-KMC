#pragma once

# include <vector>
# include <string>
# include <iostream>    // using IO functions
# include <string>      // using string
# include <fstream>
# include <sstream>

# include "file_reader.hpp"

using namespace std;

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
		virtual void simulate(int);								// perform STS KMC simulation
        void initialize_sim(int);
        virtual double get_micro_scale_sens_profile(int, int, int);

};
