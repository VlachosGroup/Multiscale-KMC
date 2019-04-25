# Multiscale-KMC

Kinetic Monte Carlo code which simulates a reaction setwork in either one or two time scales. Likelihood ratio
sensitivity analysis is used so that the sensitivties of species populations with respect to the rate constants
are computed without the need for additional runs.

![Imgur](http://i.imgur.com/5ROh9m1.png)

## Features
* Implemented in C++
* Replicate trajectories parallelized with MPI
* Input/output file interface
* Mass action kinetics, well mixed
* One and two time-scale options (two time-scale version under development)

## Getting started
* Download or clone the repository.
* Navigate to the folder with the version you want to run. Serial_run has the serial version while MPI_run has the parallel version.
* Type ```make``` to compile. This will generate an executable called MSA-KMC.x
* (optional) Type ```make clean``` to delete the object files (*.o)
* Open the input file network.in and modify it to the system you want to run. See the below for details on input file format.
* Run the executable. The parallel version must be run with ```mpiexec -n MSA-KMC.x```, where ```n``` is the number of processors.
* On squidward.che.udel.edu, it can be submitted as a batch job using ```qsub submit_serial.qs``` for the serial version or ```qsub mpi_cpp.qs``` for the parallel version.
* The code will generate output files with average species profiles and sensitivity profiles. See below for details on output file format.

## References
* [A. Hashemi, M. Nunez, P. Plechac, D.G. Vlachos, “Stochastic Averaging and Sensitivity Analysis for Two Scale Reaction Networks” Journal of Chemical Physics 144, 074104 (2016)](http://arxiv.org/abs/1509.03802)  
* [M. Núñez and D.G. Vlachos, J. Chem. Phys. 142 (4), 044108 (2015)](http://scitation.aip.org/content/aip/journal/jcp/142/4/10.1063/1.4905957)
* [Marcel Nunez Ph.D. thesis, University of Delaware, 2018](http://udspace.udel.edu/handle/19716/23604)


## Input file format

Here we tell you how to set up your KMC system in MSA_in.txt. Examples are provided in the serial and parallel run folders. Lines which begin with a # are comments. Blank spaces at the beginning and ends of lines are ignored.

| Command | Data type | Description |
| --- | --- | --- |
| Number of species | integer | Number of species in the system |
| Number of reactions | integer | Number of reactions in the system |
| Species names | string | Names of each species, used only in the output file |
| Parameter names | string | Names of each rate constant, typically numbered k1, k2, ... |
| Initial state | integer | Initial populations of each species |
| Reactions | integer | Stoichiometric matrix for the reaction network. Each line is a reaction. |
| Rate constants | double | Rate constants for each elementary step. Propensities are a based on mass action kinetics |
| Final time | double | Number of seconds (KMC time) to simulate |
| Number of trajectories | integer | Number of trajectories to use for averaging. For parallel runs, it will round this up to the nearest multiple of the number of processors. |
| write trajectory data | flag | Include this flag if you want to write output files for each individual trajectory. This will take up a lot of storage and slow down your calculation. Best used for debugging purposes. |
| two time scale | flag | Include this flag to activate two time scale mode. |
| Fast reactions | integer | Each line afterwards consists of pairs of forward and reverse reactions to be labeled as fast |


## Output file format

There are two output files.

### species_avg_out.txt

Contains a table with the population profiles for each species.

```
Time    A       B       C
0	100	0	0	
0.1	55.14	37.102	7.758
...
9.9	0.015	0.018	99.967	
10	0.009	0.023	99.968	
```

### sensitivities_out.txt

Contains a table for each species, with the sensitivity profiles for each parameter.

```
A

Time 	k1	      k2	  k3	
0	0	      0	          0	
0.1	-34.3618      -12.3738	  -3.43272
...
9.9	0.60226	      0.115044	  -0.13159	
10	-0.706654     -0.892604	  -0.05226

B

Time 	k1	      k2	  k3	
0	0	      0	          0	
0.1	28.0757	      2.15473	  -3.78631
...
9.9	-0.995014     -0.867165	  -0.156703	
10	-0.0193514    -0.0513063  -0.2276	

...	
```
