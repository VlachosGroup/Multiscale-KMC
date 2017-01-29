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
    
    // Some variables have default values. Others MUST be read from the input file

    int n_rxns;
    int n_specs;
    double t_final;                         // termination time (s)
    int N_record = 100;                     // number of time points to record
    int N_traj = 100;                       // number of replicate trajectories
    double eps;                	            // stiffness level, eps << 1
    int num_batches = 50;        	        // number of batches to use for microscale steady-state averaging in KMC_TTS
    double delta = 0.05;              	    // accuracy level of microscale averaging, delta << 1
    bool write_traj_file = false;           // set to true if the user wants data files for each individual trajectory, will take up a lot of space
    bool two_time_scale = false;            // simulation uses one time scale by default
    
    // Specify the reaction system - used for all models
    string spec_names[] = {"A","B","C"};                   // names of chemical species
    int N_0[3] = {100, 0, 0};                             // initial state
    int stoich_mat[3][3] = { {-1, 1, 0},
        {1, -1, 0},
        {0, -1, 1}};                       // stoichiometry matrix                      
    double rate_const[3] = {20.0, 30.0, 2.0};                            // rate constants of elementary reactions
    string param_names[3] = {"k1","k2","k3"};                  // names of the parameters
    
    
    int fast_rxns[2] = {1,2};                       // indices of the fast reactions
    
    
    //int n_params = n_rxns;
    
     
    
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
           
            
            if(line == "Species"){
                
                // read names of species
                
            }else if(line == "Initial state"){
                
                // read initial species populations
                
            }else if(line == "Number of species"){
                
                getline (myfile,line);
                n_specs = atoi(line.c_str());
                
            }else if(line == "Final time"){
            
                getline (myfile,line);
                t_final = atof(line.c_str());
                
            }else if(line == "Number of reactions"){
                
                getline (myfile,line);
                n_rxns = atoi(line.c_str());
                
            }else if(line == "Reactions"){
                
                // read reactions
                
            }else if(line == "Rate constants"){
                
                // read rate constants
                
            }else if(line == "write trajectory data"){      // Flag to print out data for individual trajectories
                
                write_traj_file = true;
                
            }else if(line == "two time scale"){      // Flag to turn two time scale mode on
                
                two_time_scale = true;
                
            }else if(line == "Parameter names"){
                
                // read parameter names
                
            }else if(line == "Fast reactions"){
                
                // read the fast reactions
                
            }else{                  // Unrecognized command. Throw an error.
                
                //cout << "Unrecognized command: " << line << endl;
                //return -1;
                // should also tell them the line number
                
            }
            
        }
        
        myfile.close();
    }
    
    else cout << "Unable to open file"; 
    
    cout << "There are " << n_specs << " species and " << n_rxns << " reactions." << endl;
    cout << "Simulation will take place until " << t_final << " seconds." << endl;
    
    return 0;
  
}

// For two time scale system, make sure the choice of fast reactions is physical