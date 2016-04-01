# Multiscale-KMC-with-SA
Two time scale well-mixed kinetic Monte Carlo. Multiscale likelihood ratio sensitivity analysis is done on the fly. Implemented in Matlab.

Code used in A. Hashemi, M. Nunez, P. Plechac, D.G. Vlachos, “Stochastic Averaging and Sensitivity Analysis for Two Scale Reaction Networks” Journal of Chemical Physics 144, 074104 (2016)

Includes the simulation code as well as post-processing scripts. Simualtions create data files.

To Do
- Read reaction network and rate information from an input file.
- Put modes for steady-state vs. transient

Desired functionality
- User inputs reaction network in input files
- User chooses single time-scale or two time-scale mode
- User specifies which reactions are fast and slow
- Choose simulation parameters (time horizon, sampling times, number or replicates, etc.)
- Do averageing over replicates for transient, don't do macro steady-state because that's more complicated
- Simple post-processing scripts which read the output files and output trajectories and sensitivity estimates with error estimates

- Need this to by nice and user friendly so I can quickly test different reaction networks in my thesis
- Do not include analytical solutions. Make the code general
- long term: Consider a Fortran implementation for speed
