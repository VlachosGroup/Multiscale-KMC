# MSA-KMC
Multiscale sensitivity analysis kinetic Monte Carlo

Kinetic Monte Carlo code which simulates a reaction setwork in either one or two time scales. Likelihood ratio
sensitivity analysis is used so that the sensitivties of species populations with respect to the rate constants
are computed without the need for additional runs.

### Features
* Implemented in C++
* Replicate trajectories parallelized with MPI
* Input/output file interface
* Mass action kinetics, well mixed
* One and two time-scale options (two time-scale version under development)

### Getting started
* Download or clone the repository.
* Navigate to the folder with the version you want to run. Serial_run has the serial version while MPI_run has the parallel version.
* Type ```make``` to compile. This will generate an executable called MSA-KMC.x
* (optional) Type ```make clean``` to delete the object files (*.o)
* Open the input file network.in and modify it to the system you want to run. See the [input file format](https://github.com/VlachosGroup/MSA-KMC/wiki/Input-file-format) page for details.
* Run the executable. The parallel version must be run with ```mpiexec -n MSA-KMC.x```, where ```n``` is the number of processors.
* On squidward.che.udel.edu, it can be submitted as a batch job using ```qsub submit_serial.qs``` for the serial version or ```qsub mpi_cpp.qs``` for the parallel version.
* The code will generate output files with average species profiles and sensitivity profiles. See the [output file format](https://github.com/VlachosGroup/MSA-KMC/wiki/Output-file-format) page for details.

### References
* [A. Hashemi, M. Nunez, P. Plechac, D.G. Vlachos, “Stochastic Averaging and Sensitivity Analysis for Two Scale Reaction Networks” Journal of Chemical Physics 144, 074104 (2016)](http://arxiv.org/abs/1509.03802)  
* [M. Núñez and D.G. Vlachos, J. Chem. Phys. 142 (4), 044108 (2015)](http://scitation.aip.org/content/aip/journal/jcp/142/4/10.1063/1.4905957)
* Marcel Nunez Ph.D. thesis, University of Delaware, 2017 (in preparation)

### Additional documentation
[Wiki documentation](https://github.com/VlachosGroup/MSA-KMC/wiki)