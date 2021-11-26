#include "KMC_traj.hpp"
#include <fstream>
#include <iostream> // using IO functions
#include <math.h>
#include <random>
#include <sstream>
#include <string> // using string
#include <vector>
using namespace std;

// Static variables - names of output files
string KMC_traj ::species_out_flname = "species_out.txt";
string KMC_traj ::traj_deriv_out_flname = "traj_deriv_out.txt";

KMC_traj ::KMC_traj()
    : in_data(), writer_spec(), writer_SA() {} // empty constructor

void KMC_traj ::simulate(int rand_seed) { // Execute simulation

  srand(rand_seed); // Set the random seed

  // Some additional non-class varaiables
  vector<double> prop_cum;
  prop_cum.resize(in_data.n_rxns);
  double asum;
  double r_rxn_choose; // random number between 0 and 1 used to choose which
                       // reaction to fire
  double
      r_timestep; // random mumber between 0 and 1 used to choose the time step

  initialize_sim(rand_seed);

  // Start KMC loop
  while (t < in_data.t_final) {

    /*
    ============================ Compute quantities for current step
    ============================
    */

    // Generate random numbers
    r_rxn_choose = ((double)rand() / (RAND_MAX));
    r_timestep = ((double)rand() / (RAND_MAX));

    // Compute reaction propensities
    for (int i = 0; i < in_data.n_rxns; i++) {

      // Use mass action kinetics equation, a = k * [A]^ma * [B]^mb * ...
      props[i] = in_data.rate_const[i];
      for (int j = 0; j < in_data.n_specs; j++) {
        if (in_data.stoich_mat[i][j] <
            0) { // if this species is a reactant, use it to compute the rate
          props[i] = props[i] * pow(N[j], -in_data.stoich_mat[i][j]);
        }
      }

      // Fill in derivatives for each parameter
      for (int k = 0; k < in_data.n_params; k++) {
        if (i == k) {

          prop_ders[i][k] = 1.0;

          for (int j = 0; j < in_data.n_specs; j++) {
            if (in_data.stoich_mat[i][j] <
                0) { // if this species is a reactant, use it to compute the
                     // rate
              prop_ders[i][k] =
                  prop_ders[i][k] * pow(N[j], -in_data.stoich_mat[i][j]);
            }
          }

        } else {
          prop_ders[i][k] = 0;
        }
      }
    }

    // Add derivatives for each parameter
    for (int k = 0; k < in_data.n_params; k++) {
      prop_ders_sum[k] = 0;
      for (int i = 0; i < in_data.n_rxns; i++) {
        prop_ders_sum[k] += prop_ders[i][k];
      }
    }

    // Sum the propensities
    asum = 0;
    for (int i = 0; i < in_data.n_rxns; i++) {
      asum += props[i];
    }

    // If all propensities are 0, then exit the while loop
    if (asum == 0) {
      break;
    }

    // Choose which reaction to fire
    for (int i = 0; i < in_data.n_rxns; i++) {
      prop_cum[i] = 0;
      for (int j = 0; j <= i; j++) {
        prop_cum[i] += props[j] / asum;
      }
    }

    rxn_to_fire_ind = 0;
    if (r_rxn_choose == 1) { // If the random number ends up being zero, choose
                             // the last reaction with nonzero propensity
      for (int i = 0; i < in_data.n_rxns; i++) {
        if (props[i] > 0) {
          rxn_to_fire_ind = i;
        }
      }
    } else {
      while (prop_cum[rxn_to_fire_ind] < r_rxn_choose) {
        rxn_to_fire_ind += 1;
      }
    }

    if (r_rxn_choose == 1) {
      cout << rxn_to_fire_ind << endl;
    }

    // Compute time step
    dt = log(1 / r_timestep) / asum;
    if (not std::isfinite(dt)) {
      break;
    }

    // Record the current state as long as time >= t_sample
    while (t >= in_data.t_rec[ind_rec]) {
      record_stats();
    }

    /*
    ============ Fire reaction and update system ==============
    */

    // Update species populations
    for (int j = 0; j < in_data.n_specs; j++) {
      N[j] += in_data.stoich_mat[rxn_to_fire_ind][j];
    }

    // Update clock
    t_prev = t;
    t = t + dt;

    // Update trajectory derivatives
    for (int k = 0; k < in_data.n_params; k++) {
      if (props[rxn_to_fire_ind] == 0) {
        cout << "Propensity is zero" << endl;
        throw 20;
      }
      W[k] += prop_ders[rxn_to_fire_ind][k] /
              props[rxn_to_fire_ind]; // contribution from reaction that fires
      W[k] -= prop_ders_sum[k] * dt;  // contribution from time step
      if (not std::isfinite(W[k])) {
        cout << "W is NaN" << endl;
        cout << "The chosen reaction is " << endl;
        cout << rxn_to_fire_ind << endl;
        cout << "The propensity is " << endl;
        cout << props[rxn_to_fire_ind] << endl;
        cout << "The sum of propensity derivatives is" << endl;
        cout << prop_ders_sum[k] << endl;
        cout << "The time step is " << endl;
        cout << dt << endl;
        cout << "Random number for choosing reaction" << endl;
        cout << r_rxn_choose << endl;
        cout << "Random number for time step" << endl;
        cout << r_timestep << endl;
        throw 20;
      }
    }
  }

  // Fill in recording times that were missed
  while (ind_rec < in_data.N_record) {
    record_stats();
  }

  // Close output files
  if (in_data.write_traj_files) {
    writer_spec.close();
    writer_SA.close();
  }
}

void KMC_traj ::initialize_sim(int rand_seed) {

  /*
  ============================ Initialize variables ============================
  */

  t = 0;       // KMC clock
  t_prev = 0;  // KMC time of previous step
  ind_rec = 0; // time point

  // Initialize species populations
  N.resize(in_data.n_specs);
  for (int k = 0; k < in_data.n_specs; k++) {
    N[k] = in_data.N_0[k];
  }

  props.resize(in_data.n_rxns);

  // Derivative of each propensity with respect to each parameter
  prop_ders_sum.resize(
      in_data.n_rxns); // derivative of the sum of all propensities with respect
                       // to each parameter
  prop_ders.resize(in_data.n_rxns); // derivative of the sum of all propensities
                                    // with respect to each parameter
  for (int i = 0; i < in_data.n_rxns; i++) {
    prop_ders[i].resize(in_data.n_params);
  }

  // Initialize trajectory derivatives
  W.resize(in_data.n_params); // trajectory derivatives
  for (int k = 0; k < in_data.n_params; k++) {
    W[k] = 0;
  }

  // Vector for recording species profiles
  spec_profile.resize(in_data.N_record);
  for (int i = 0; i < in_data.N_record; i++) {
    spec_profile[i].resize(in_data.n_specs);
  }

  // Vector for recording trajectory derivatives
  traj_deriv_profile.resize(in_data.N_record);
  for (int i = 0; i < in_data.N_record; i++) {
    traj_deriv_profile[i].resize(in_data.n_params);
  }

  // Initialize all entries as zero
  for (int i = 0; i < in_data.N_record; i++) {
    for (int j = 0; j < in_data.n_params; j++) {
      traj_deriv_profile[i][j] = 0;
    }
  }

  /*
  ============== Open files to record data for the trajectory ==============
  */

  if (in_data.write_traj_files) {

    // Create a folder for all trajectory data files
    // Name the file specific to this trajectory. Put random seed in the file
    // name.

    ostringstream fname1;
    fname1 << KMC_traj::species_out_flname << "_" << rand_seed << ".out";

    writer_spec.open(fname1.str());
    if (!writer_spec) {
      cout << "Error opening file" << endl;
    }

    // Write header
    writer_spec << "Time \t";
    for (int j = 0; j < in_data.n_specs; j++) {
      writer_spec << in_data.spec_names[j] << "\t";
    }
    writer_spec << endl;

    ostringstream fname2;
    fname2 << KMC_traj::traj_deriv_out_flname << "_" << rand_seed << ".out";

    writer_SA.open(fname2.str());
    if (!writer_SA) {
      cout << "Error opening file" << endl;
    }

    // Write header
    writer_SA << "Time \t";
    for (int i = 0; i < in_data.n_params; i++) {
      writer_SA << in_data.param_names[i] << "\t";
    }
    writer_SA << endl;
  }
}

/*
========= Record the current state =========
*/
void KMC_traj ::record_stats() {

  t_trunc = in_data.t_rec[ind_rec] - t_prev;

  // Record data in output files
  if (in_data.write_traj_files) {

    // Record the species numbers
    writer_spec << in_data.t_rec[ind_rec] << "\t";
    for (int j = 0; j < in_data.n_specs; j++) {
      writer_spec << N[j] << "\t";
    }
    writer_spec << endl;

    // Record trajectory derivatives
    writer_SA << in_data.t_rec[ind_rec] << "\t";
    for (int k = 0; k < in_data.n_params; k++) {
      writer_SA << W[k] - prop_ders_sum[k] * t_trunc << "\t";
    }
    writer_SA << endl;
  }

  // Record data in variables

  // Record species populations
  for (int j = 0; j < in_data.n_specs; j++) {
    spec_profile[ind_rec][j] = (double)N[j];
  }

  // Record trajectory derivatives
  for (int k = 0; k < in_data.n_params; k++) {
    traj_deriv_profile[ind_rec][k] = W[k] - prop_ders_sum[k] * t_trunc;
  }

  ind_rec += 1;
}

double KMC_traj ::get_micro_scale_sens_profile(int i, int j, int k) {
  return 0;
}
