#pragma once

#include <fstream>
#include <iostream> // using IO functions
#include <sstream>
#include <string>
#include <string> // using string
#include <vector>

#include "KMC_traj.hpp"
#include "KMC_traj_TTS.hpp"
#include "file_reader.hpp"

using namespace std;

/*!
Class to collect statistics from multiple KMC trajectories
*/
class Traj_stats {

private:
  static string species_avgs_out_flname; // "spec_avgs.out"
  static string SA_out_flname;           // "SA.out"

  vector<vector<double>> spec_profiles_averages; // species averages
  vector<vector<double>> traj_deriv_avgs; // trajectory derivative averages
  vector<vector<vector<double>>> sensitivities; // sensities
  vector<vector<vector<double>>> microscale_sensitivity_contributions;

public:
  file_reader in_data; // input data from the input file

  Traj_stats(); // Empty class constructor
  void initialize_stats();
  void run_simulations();
  void finalize_stats();
  void write_spec_avg_output();
  void write_sensitivity_output();
};
