# include <iostream>    // using IO functions
# include <string>      // using string
# include <vector>
# include <fstream>
# include <random>
# include <math.h>
# include <algorithm> // for std::find
# include <iterator> // for std::begin, std::end

# include <typeinfo>
# include <sstream>
# include <stdio.h>      /* printf, fgets */
# include <stdlib.h>     /* atoi */
# include "file_reader.hpp"
using namespace std;


// Helper function which removes white spaces from the beginning and end of a string
string trim_string(const string& string_in){       // remove white spaces from the beginning and the end

    int strBegin = string_in.find_first_not_of(" ");
    if (strBegin == std::string::npos)
        return ""; // no content

    int strEnd = string_in.find_last_not_of(" ");
    int strRange = strEnd - strBegin + 1;

    return string_in.substr(strBegin, strRange);
}


file_reader :: file_reader(){}      // empty constructor

// Constructor for the file_reader class
file_reader :: file_reader(string flname){


	ifstream myfile(flname);
	string line;

	if (myfile.is_open())
	{

		while ( getline (myfile,line) )
		{


			line = trim_string(line);       // Trim empty spaces from the beginning and end

			if(line == "" or line[0] == '#')    // skip empty lines and comments
			{
				continue;
			}


			if(line == "Number of species"){

				getline (myfile,line);
				n_specs = atoi(line.c_str());

                N_0.resize(n_specs);                       // allocate memory for the initial state

			}else if(line == "Number of reactions"){

				getline (myfile,line);
				n_rxns = atoi(line.c_str());
				n_params = n_rxns;

				rate_const.resize(n_rxns);

            }else if(line == "random seed"){

				getline (myfile,line);
				rand_seed = atoi(line.c_str());

			}else if(line == "Species names"){

				for(int spec_ind = 0; spec_ind < n_specs; spec_ind++){
					getline (myfile,line);
					line = trim_string(line);
                    spec_names.push_back(line);
				}

			}else if(line == "Initial state"){

				for(int spec_ind = 0; spec_ind < n_specs; spec_ind++){
					getline (myfile,line);
					line = trim_string(line);
					N_0[spec_ind] = atoi(line.c_str());
				}

			}else if(line == "Final time"){

				getline (myfile,line);
				line = trim_string(line);
				t_final = atof(line.c_str());

			}else if(line == "Reactions"){          // read in stoichiometric matrix


				for(int rxn_ind = 0; rxn_ind < n_rxns; rxn_ind++){

					getline (myfile,line);

					string delimiter = " ";
					size_t pos = 0;
					string token;

                    vector<int> rxn_stoich;

					for(int spec_ind = 0; spec_ind < n_specs; spec_ind++){
						line = trim_string(line);
						pos = line.find(delimiter);
						token = line.substr(0, pos);
                        rxn_stoich.push_back(atoi(token.c_str()));
						line.erase(0, pos + delimiter.length());
					}

                    stoich_mat.push_back(rxn_stoich);
				}

			}else if(line == "Rate constants"){

				for(int rxn_ind = 0; rxn_ind < n_rxns; rxn_ind++){
					getline (myfile,line);
					line = trim_string(line);
					rate_const[rxn_ind] = atof(line.c_str());
				}

			}else if(line == "write trajectory data"){      // Flag to print out data for individual trajectories

				write_traj_files = true;

			}else if(line == "two time scale"){      // Flag to turn two time scale mode on

				two_time_scale = true;

			}else if(line == "Number of fast reaction pairs"){      // Flag to turn two time scale mode on

				getline (myfile,line);
				n_fast_pairs = atoi(line.c_str());

                n_fast_rxns = n_fast_pairs * 2;
                n_slow_rxns = n_rxns - n_fast_rxns;
                fast_rxns.resize(n_fast_rxns);
                slow_rxns.resize(n_slow_rxns);
                fast_pairs.resize(n_fast_rxns);
                for(int i = 0; i < n_fast_rxns; i++){
                    fast_pairs[i].resize(2);
                }

            }else if(line == "Fast reaction pairs"){      // Flag to turn two time scale mode on

                int fwd_rxn_ind;
                int rev_rxn_ind;

				for(int pair_ind = 0; pair_ind < n_fast_pairs; pair_ind++){

					getline (myfile,line);

					string delimiter = " ";
					size_t pos = 0;
					string token;

                    // Read forward reaction
					line = trim_string(line);
					pos = line.find(delimiter);
					token = line.substr(0, pos);
                    fwd_rxn_ind = atoi(token.c_str()) - 1;
                    fast_pairs[2 * pair_ind][0] = fwd_rxn_ind;
                    fast_pairs[2 * pair_ind + 1][1] = fwd_rxn_ind;
					line.erase(0, pos + delimiter.length());

                    // Read reverse reaction
                    line = trim_string(line);
					pos = line.find(delimiter);
					token = line.substr(0, pos);
                    rev_rxn_ind = atoi(token.c_str()) - 1;
                    fast_pairs[2 * pair_ind][1] = rev_rxn_ind;
                    fast_pairs[2 * pair_ind + 1][0] = rev_rxn_ind;
					line.erase(0, pos + delimiter.length());

                    fast_rxns[2 * pair_ind] = fwd_rxn_ind;
                    fast_rxns[2 * pair_ind + 1] = rev_rxn_ind;

				}

                // Populate list of slow reactions
                int slow_ind = 0;
                for(int i = 0; i < n_rxns; i++){

                    if ( find(begin(fast_rxns), end(fast_rxns), i) == end(fast_rxns) ){     // If it is not a fast reaction, it is a slow reaction
                        slow_rxns[slow_ind] = i;
                        slow_ind ++;
                    }

                }


            }else if(line == "Number of microscale averaging steps"){

				getline (myfile,line);
				n_micro_steps = atoi(line.c_str());

            }else if(line == "Parameter names"){

				for(int param_ind = 0; param_ind < n_params; param_ind++){
					getline (myfile,line);
					line = trim_string(line);
                    param_names.push_back(line);
				}

            }else if(line == "Number of trajectories"){

				getline (myfile,line);
				N_traj = atoi(line.c_str());

            }else if(line == "Number of time points"){

				getline (myfile,line);
				N_record = atoi(line.c_str()) + 1;

			}else if(line == "Fast reactions"){

				// read the fast reactions

			}else{                  // Unrecognized command. Throw an error. Should also tell them the line number

				cout << "Unrecognized command: " << line << endl;
			}

		}
	}else{
		cout << "Unable to open file" << endl;
	}

	myfile.close();


    // Fill in recording times
    t_rec.resize(N_record);
	for(int i = 0; i < N_record; i++){
		t_rec[i] = (double) i / (N_record-1) * t_final;
	}


}
