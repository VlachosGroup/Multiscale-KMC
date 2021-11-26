#include "Traj_stats.hpp"
#include "file_reader.hpp"

int main() {

  cout.precision(6);
  file_reader fr("network.in"); // Read input file
  Traj_stats sim;
  sim.in_data = fr;      // assign input
  sim.run_simulations(); // run KMC simulations, will execute differently
                         // depending on fr.two_time_scale

  return 0;
}
