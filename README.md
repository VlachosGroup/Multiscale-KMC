# Multiscale-KMC-with-SA
Two time scale well-mixed kinetic Monte Carlo. Multiscale likelihood ratio sensitivity analysis is done on the fly. Implemented in Matlab.

Journal Articles
- Code used in A. Hashemi, M. Nunez, P. Plechac, D.G. Vlachos, “Stochastic Averaging and Sensitivity Analysis for Two Scale Reaction Networks” Journal of Chemical Physics 144, 074104 (2016)
- Two-time scale sensitivity analysis method developed in M. Núñez and D.G. Vlachos, J. Chem. Phys. 142 (4), 044108 (2015)

Includes the simulation code as well as post-processing scripts. Simualtion codes create data files which are read and analyzed by the post-processing scripts. The ODE scripts combine both into a single m files because the ODEs are fast.

To Do
- Read reaction network and rate information from an input file.
- Put modes for steady-state vs. transient and STS vs. TTS
- Include a latex document with a simple example worked out (with data). There is an example in Research_Files/Thesis/MultiscaleSA
- Complete automaton for arbitrary networks in ODE_TTS
- Make the graph analysis in each function consistent, complete, and general for arbitrary networks. 
	May want to just have each function output a trajectory of population averages and sensitivities.
	
Desired functionality
- Input reaction network in input files
- Choose between deterministic (ODE) and stochastic (KMC)
- Choose single time-scale (STS) or two time-scale mode (TTS)
- Specify which reactions are fast and slow
- Choose simulation parameters (time horizon, sampling times, number or replicates, etc.)
- Averageing over replicates for transient (don't do macro steady-state because that's more complicated)
- Simple post-processing scripts which read the output files and output trajectories and sensitivity estimates with error estimates

Future Work
- Port to a Fortran implementation for speed
