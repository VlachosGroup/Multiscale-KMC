# include <iostream>    // using IO functions
# include <string>      // using string
# include <vector>
# include <fstream>
# include <random>
# include <math.h>

# include <typeinfo>
# include <sstream>
# include <stdio.h>      /* printf, fgets */
# include <stdlib.h>     /* atoi */
using namespace std;

string trim_string(const string& string_in){       // remove white spaces from the beginning and the end
    
    int strBegin = string_in.find_first_not_of(" ");
    if (strBegin == std::string::npos)
        return ""; // no content

    int strEnd = string_in.find_last_not_of(" ");
    int strRange = strEnd - strBegin + 1;

    return string_in.substr(strBegin, strRange);
}


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


		file_reader(string flname){
			
			// Set default values                      
			N_record = 100;                    
			N_traj = 1000;
            num_batches = 50;        	       
			delta = 0.05;                           	   
			
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
						
                        N_0.reserve(n_specs);                       // allocate memory for the initial state          
					
					}else if(line == "Number of reactions"){
						
						getline (myfile,line);
						n_rxns = atoi(line.c_str());
						n_params = n_rxns;
						
						rate_const.reserve(n_rxns);
					
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
						N_record = atoi(line.c_str());
                    
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

		}
        
		if(two_time_scale){
            // prepare the TTS variables: fast and slow stoichiometry matrices, n_rxns_fast and n_rxns_slow
        }
        
	};

// Test driver function
int main() {
	
	file_reader fr("network.in");            // Read input file
	
    // If it is STS, make a Traj_stats_STS object and run it
    
    // If it is TTS, make a Traj_stats_TTS object and run and process it slightly differently

	return 0;
}