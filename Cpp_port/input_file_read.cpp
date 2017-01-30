// Read input.txt

# include <iostream>
# include <vector>
# include <fstream>
# include <string>
# include <typeinfo>
# include <sstream>
# include <stdio.h>      /* printf, fgets */
# include <stdlib.h>     /* atoi */
using namespace std;

int main () {
    
    /*
    =========== Declare variables to be filled in by reading the input file 
    */
    
    // Variables required for KMC simulation
    int n_rxns;
    int n_specs;
    int n_params;
    double t_final;                         // termination time (s)
    int N_record = 100;                     // number of time points to record
    int N_traj = 1000;                       // number of replicate trajectories
    double eps;                	            // stiffness level, eps << 1
    bool write_traj_files = false;           // flag to write data files for each individual trajectory, will take up a lot of space
    string* spec_names = NULL;              // names of chemical species
    int* N_0 = NULL;                        // initial state
    int** stoich_mat = NULL;                // stoichiometry matrix                      
    double* rate_const = NULL;              // rate constants of elementary reactions
    string* param_names = NULL;             // names of the parameters
    
    // Two time scale variables
    bool two_time_scale = false;            // simulation uses one time scale by default
    int* fast_rxns = NULL;                  // indices of the fast reactions
    int num_batches = 50;        	        // number of batches to use for microscale steady-state averaging in KMC_TTS
    double delta = 0.05;              	    // accuracy level of microscale averaging, delta << 1
     
    
    /*
    ====================== Read the input file ======================
    */
    
    string line;
    
    ifstream myfile ("input.txt");
  
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
            
            if(line.find_first_not_of(' ') == std::string::npos or line[0] == '#')    // skip empty lines and comments
            {
                continue;
            }
           
            
            if(line == "Number of species"){
                
                getline (myfile,line);
                n_specs = atoi(line.c_str());
                
                spec_names = new string[n_specs];       // allocate memory for the species names
                N_0 = new int[n_specs];                 // allocate memory for the initial state
            
            }else if(line == "Number of reactions"){
                
                getline (myfile,line);
                n_rxns = atoi(line.c_str());
                n_params = n_rxns;
                
                param_names = new string[n_params];       // allocate memory for the species names
                rate_const = new double[n_rxns];
            
            }else if(line == "Species names"){
                
                for(int spec_ind = 0; spec_ind < n_specs; spec_ind++){
                    getline (myfile,line);
                    spec_names[spec_ind] = line;
                }
                
            }else if(line == "Initial state"){
                
                for(int spec_ind = 0; spec_ind < n_specs; spec_ind++){
                    getline (myfile,line);
                    N_0[spec_ind] = atoi(line.c_str());
                }
                 
            }else if(line == "Final time"){
            
                getline (myfile,line);
                t_final = atof(line.c_str());
                
            }else if(line == "Reactions"){          // read in stoichiometric matrix
                
                stoich_mat = new int*[n_rxns];          // allocate stoichiometric matrix
                for(int i = 0; i < n_rxns; ++i){
                    stoich_mat[i] = new int[n_specs];
                }
                    
                for(int rxn_ind = 0; rxn_ind < n_rxns; rxn_ind++){
                    
                    getline (myfile,line);
                    string delimiter = " ";
                    //string token = line.substr(0, line.find(" ")); // token is "scott"
                    size_t pos = 0;
                    string token;
                    
                    for(int spec_ind = 0; spec_ind < n_specs; spec_ind++){
                        pos = line.find(delimiter);
                        token = line.substr(0, pos);
                        stoich_mat[rxn_ind][spec_ind] = atoi(token.c_str());
                        line.erase(0, pos + delimiter.length());
                    } 
                    
                }
                
            }else if(line == "Rate constants"){
                
                for(int rxn_ind = 0; rxn_ind < n_rxns; rxn_ind++){
                    getline (myfile,line);
                    rate_const[rxn_ind] = atof(line.c_str());
                }
                
            }else if(line == "write trajectory data"){      // Flag to print out data for individual trajectories
                
                write_traj_files = true;
                
            }else if(line == "two time scale"){      // Flag to turn two time scale mode on
                
                two_time_scale = true;
                
            }else if(line == "Parameter names"){
                
                for(int param_ind = 0; param_ind < n_params; param_ind++){
                    getline (myfile,line);
                    param_names[param_ind] = line;
                }
                
            }else if(line == "Fast reactions"){
                
                // read the fast reactions
                
            }else{                  // Unrecognized command. Throw an error. Should also tell them the line number
                
                cout << "Unrecognized command: " << line << endl;
                return -1;
                
            }
            
        }
        
        myfile.close();
    }
    
    else cout << "Unable to open file"; 
    
    cout << "Successfully read input file" << endl;
    
    return 0;
  
}