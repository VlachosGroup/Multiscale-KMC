# include <iostream>    // using IO functions
# include <string>      // using string
# include <vector>
# include <fstream>
# include <random>
# include <math.h>
# include "myHeader.h"
#ifdef MPI
# include <mpi.h>
#endif
using namespace std;

string Traj_stats :: species_avgs_out_flname = "species_avgs_out.txt";
string Traj_stats :: SA_out_flname = "sensitivities_out.txt";

Traj_stats :: Traj_stats() : in_data(){}      // Empty class constructor


void Traj_stats :: initialize_stats(){

    // Resize vectors to be able to hold the necessary data

    spec_profiles_averages.resize(in_data.N_record);
    for (int i = 0; i < in_data.N_record; i++){
        spec_profiles_averages[i].resize(in_data.n_specs);}

    traj_deriv_avgs.resize(in_data.N_record);
    for (int i = 0; i < in_data.N_record; i++){
        traj_deriv_avgs[i].resize(in_data.n_params);}

    sensitivities.resize(in_data.N_record);
    for (int i = 0; i < in_data.N_record; i++){

        sensitivities[i].resize(in_data.n_specs);
        for(int j = 0; j < in_data.n_specs; j++){

            sensitivities[i][j].resize(in_data.n_params);
        }
    }


    // Intialize statistics

    for (int i = 0; i < in_data.N_record; ++i){
        for(int j = 0; j < in_data.n_specs; j++){
            spec_profiles_averages[i][j] = 0;
        }
    }

    for (int i = 0; i < in_data.N_record; ++i){
        for(int j = 0; j < in_data.n_params; j++){
            traj_deriv_avgs[i][j] = 0;
        }
    }

    for (int i = 0; i < in_data.N_record; ++i){
        for(int j = 0; j < in_data.n_specs; j ++){
            for(int k = 0; k < in_data.n_params; k++){
                sensitivities[i][j][k] = 0;
            }
        }
    }

}


void Traj_stats :: finalize_stats(){

    /*
    ============ Finalize statistical analysis ==============
    */

    // Normalize statistics by the number of trajectories
    for (int i = 0; i < in_data.N_record; ++i){

        for(int j = 0; j < in_data.n_specs; j++){
            spec_profiles_averages[i][j] = spec_profiles_averages[i][j] / in_data.N_traj;

            for(int k = 0; k < in_data.n_params; k++){
                sensitivities[i][j][k] = sensitivities[i][j][k] / in_data.N_traj;

                if ( std::isnan( sensitivities[i][j][k] ) ){
                    cout << "Detected a NaNs" << endl;
                    throw 20;
                }
            }
        }

        for(int j = 0; j < in_data.n_params; j++){
            traj_deriv_avgs[i][j] = traj_deriv_avgs[i][j] / in_data.N_traj;
        }
    }


    // Subtract 2nd term of the covariance, renormalize by N_traj-1
    for (int i = 0; i < in_data.N_record; ++i){

        for(int j = 0; j < in_data.n_specs; j++){

            for(int k = 0; k < in_data.n_params; k++){

                //cout << endl;
                //cout << sensitivities[i][j][k] << endl;
                sensitivities[i][j][k] -= spec_profiles_averages[i][j] * traj_deriv_avgs[i][k];
                sensitivities[i][j][k] = sensitivities[i][j][k] * in_data.N_traj / (in_data.N_traj - 1);
                //cout << sensitivities[i][j][k] << endl;
            }
        }

    }

}


/*
============================ Run simulations to gather data ============================
*/

#ifndef MPI

void Traj_stats :: run_simulations(){

    initialize_stats();

    for(int traj_ind = 0; traj_ind < in_data.N_traj; traj_ind++){

        // Create and run a KMC simulation

        //cout << traj_ind + 1 << " / " << in_data.N_traj << endl;              // Print the trajectory number

        KMC_traj* run = NULL;  // initalize the pointer

        if(in_data.two_time_scale){     // Two time scale
            run = new KMC_traj_TTS;
        }else{                          // Single time scale
            run = new KMC_traj;
        }

        run->in_data = in_data;             // Copy input file data to the trajectory object
        run->simulate(in_data.rand_seed + traj_ind);

        // Add to statistical running counts
        for (int i = 0; i < in_data.N_record; ++i){

            for(int j = 0; j < in_data.n_specs; j++){
                spec_profiles_averages[i][j] += run->spec_profile[i][j];

                for(int k = 0; k < in_data.n_params; k++){  // For TTS, add extra contribution from microscale averaging
                    sensitivities[i][j][k] += run->spec_profile[i][j] * run->traj_deriv_profile[i][k] + run->get_micro_scale_sens_profile(i, j, k);
                        if ( not std::isfinite( sensitivities[i][j][k] ) ){
                            cout << traj_ind << " has NaNs" << endl;
                            cout << in_data.rand_seed + traj_ind << " is the random seed." << endl;
                        }
                }
            }

            for(int j = 0; j < in_data.n_params; j++){
                traj_deriv_avgs[i][j] += run->traj_deriv_profile[i][j];
            }
        }

    }


    finalize_stats();

    // Write output files
    write_spec_avg_output();
    write_sensitivity_output();

    cout << "Simulation complete." << endl;
}


#else


void Traj_stats :: run_simulations(){       // add some if statements which will change bahavior depending on if it is STS or TTS

    /*
    ============================ Run simulations to gather data ============================
    */

    initialize_stats();

    // Set up MPI variables

    int id;
    int ierr;
    int p;
    int Npp;
    //int argc;
    //char *argv[];

    ierr = MPI_Init ( NULL, NULL );                 //  Initialize MPI.
    ierr = MPI_Comm_size ( MPI_COMM_WORLD, &p );      //  Get the number of processes.
    ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id );     //  Get the individual process ID.



    Npp = in_data.N_traj / p + 1;
    in_data.N_traj = Npp * p;       // round up to the nearest multiple of the number of processors

    if(id==0){
        cout << in_data.N_traj << " replicate trajectories will be used." << endl;
    }


    // Make a data array for each processor to be combined later

    int size1 = in_data.N_record * in_data.n_specs;
    int size2 = in_data.N_record * in_data.n_params;
    int size3 = in_data.N_record * in_data.n_specs * in_data.n_params;

    double spec_profiles_averages_arr[size1];
    double traj_deriv_avgs_arr[size2];
    double sensitivities_arr[size3];

    for(int i = 0; i < size1; i++){
        spec_profiles_averages_arr[i] = 0;
    }

    for(int i = 0; i < size2; i++){
        traj_deriv_avgs_arr[i] = 0;
    }

    for(int i = 0; i < size3; i++){
        sensitivities_arr[i] = 0;
    }

    int ind1;
    int ind2;
    int ind3;

    /*
    ============ run KMC simulations to gather data ==============
    */

    for(int traj_ind = 0; traj_ind < Npp; traj_ind++){

        KMC_traj* run = NULL;  // initalize the pointer

        if(in_data.two_time_scale){     // Two time scale
            run = new KMC_traj_TTS;
        }else{                          // Single time scale
            run = new KMC_traj;
        }


        run->in_data = in_data;             // Copy input file data to the trajectory object
        run->simulate(in_data.rand_seed + traj_ind + id * Npp);

        // Add to statistical running counts

        ind1 = 0;
        ind2 = 0;
        ind3 = 0;

        for (int i = 0; i < in_data.N_record; ++i){

            for(int j = 0; j < in_data.n_specs; j++){
                spec_profiles_averages_arr[ind1] += run->spec_profile[i][j];
                ind1++;

                for(int k = 0; k < in_data.n_params; k++){  // For TTS, add extra contribution from microscale averaging
                    sensitivities_arr[ind3] += run->spec_profile[i][j] * run->traj_deriv_profile[i][k] + run->get_micro_scale_sens_profile(i, j, k);
                    ind3++;
                }
            }

            for(int j = 0; j < in_data.n_params; j++){
                traj_deriv_avgs_arr[ind2] += run->traj_deriv_profile[i][j];
                ind2++;
            }
        }

    }

    /*
    ============ Collect data from the replicate simulations ==============
    */

    // For final analysis, gather data into the big arrays
    double spec_profiles_averages_big_arr[p * size1];
    double traj_deriv_avgs_big_arr[p * size2];
    double sensitivities_big_arr[p * size3];

    MPI_Gather(spec_profiles_averages_arr, size1, MPI_DOUBLE, spec_profiles_averages_big_arr, size1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(traj_deriv_avgs_arr, size2, MPI_DOUBLE, traj_deriv_avgs_big_arr, size2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(sensitivities_arr, size3, MPI_DOUBLE, sensitivities_big_arr, size3, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    // Reshape the data into multidimensional arrays on processor 0
    if(id==0){

        ind1 = 0;
        ind2 = 0;
        ind3 = 0;

        //for

        // Normalize statistics by the number of trajectories

        for(int proc_id = 0; proc_id < p; proc_id++){

            for (int i = 0; i < in_data.N_record; ++i){

                for(int j = 0; j < in_data.n_specs; j++){

                    spec_profiles_averages[i][j] += spec_profiles_averages_big_arr[ind1];
                    ind1++;

                    for(int k = 0; k < in_data.n_params; k++){

                        sensitivities[i][j][k] += sensitivities_big_arr[ind3];
                        ind3++;

                    }
                }

                for(int j = 0; j < in_data.n_params; j++){
                        traj_deriv_avgs[i][j] += traj_deriv_avgs_big_arr[ind2];
                        ind2++;
                }


            }
        }


        finalize_stats();
        write_spec_avg_output();                // Write output files
        write_sensitivity_output();
        cout << "Simulation completed successfully" << endl;
    }

    MPI_Finalize ( );     //  Terminate MPI.

}


#endif


/*
=========== Print species population averages into an output file ===========
*/
void Traj_stats :: write_spec_avg_output(){

    ofstream writer_spec_avgs;

    writer_spec_avgs.open(Traj_stats :: species_avgs_out_flname);
    if(! writer_spec_avgs){
        cout << "Error opening file" << endl;
    }

    // Write header
    writer_spec_avgs << "Species mean populations versus time" << endl;
    writer_spec_avgs << endl;
    writer_spec_avgs << "Time \t";
    for(int j = 0; j < in_data.n_specs; j++){
        writer_spec_avgs << in_data.spec_names[j] << "\t";
    }
    writer_spec_avgs << endl;

    for(int time_ind = 0; time_ind < in_data.N_record; time_ind ++){

        writer_spec_avgs << in_data.t_rec[time_ind] << "\t";

        for(int spec_ind = 0; spec_ind < in_data.n_specs; spec_ind++){
            writer_spec_avgs << spec_profiles_averages[time_ind][spec_ind] << "\t"; // print the mean population of this species at this time
        }

        writer_spec_avgs << endl;
    }

    writer_spec_avgs.close();

}


/*
===========Print sensitivities into an output file ===========
*/
void Traj_stats :: write_sensitivity_output(){

    ofstream writer_sensitivities;

    writer_sensitivities.open(Traj_stats :: SA_out_flname);
    if(! writer_sensitivities){
        cout << "Error opening file" << endl;
    }

    // Write header
    writer_sensitivities << "Sensitivities for each species and rate constant versus time" << endl << endl;

    // Loop over species
    for(int i=0; i < in_data.n_specs; i++){

        writer_sensitivities << in_data.spec_names[i] << endl << endl;

        // Write header for this species
        writer_sensitivities << "Time \t";
        for(int k = 0; k < in_data.n_params; k++){
            writer_sensitivities << in_data.param_names[k] << "\t";
        }
        writer_sensitivities << endl;

        // Loop over parameters to print out
        for(int time_ind = 0; time_ind < in_data.N_record; time_ind ++){

            writer_sensitivities << in_data.t_rec[time_ind] << "\t";

            for(int param_ind = 0; param_ind < in_data.n_params; param_ind++){
                writer_sensitivities << sensitivities[time_ind][i][param_ind] * in_data.rate_const[param_ind] << "\t";      // normalize by rate constant value
            }

            writer_sensitivities << endl;
        }

        // Put a barrier in between species data
        if(i < in_data.n_specs - 1){
            writer_sensitivities << endl;
            writer_sensitivities << "=======================================" << endl;
            writer_sensitivities << endl;
        }

    }

    writer_sensitivities.close();

}
