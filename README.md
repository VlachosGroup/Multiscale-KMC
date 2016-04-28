# Multiscale-KMC-with-SA
Two time-scale well-mixed kinetic Monte Carlo. Multiscale likelihood ratio sensitivity analysis is done on the fly. Implemented in Matlab. Single time-scale and ODE analogues are also included.

Journal Articles
- Code used in A. Hashemi, M. Nunez, P. Plechac, D.G. Vlachos, “Stochastic Averaging and Sensitivity Analysis for Two Scale Reaction Networks” Journal of Chemical Physics 144, 074104 (2016)
- Two-time scale sensitivity analysis method developed in M. Núñez and D.G. Vlachos, J. Chem. Phys. 142 (4), 044108 (2015)

Please see wiki page for workflow instructions.

Workflow
- Input reaction network in input files
- Specify which reactions are fast and slow
- Choose numerical parameters (time horizon, sampling times, number or replicates, etc.)
- Choose a model: deterministic (ODE) vs. stochastic (KMC), single time-scale (STS) vs. two time-scale mode (TTS)
- Averageing over replicates for transient (don't do macro steady-state because that's more complicated)
- Simple post-processing scripts which read the output files and output trajectories and sensitivity estimates with error estimates

To Do
- Read reaction network and rate information from an input file.
- Include a latex document with a simple example worked out (with data). There is an example in Research_Files/Thesis/MultiscaleSA
- Make the graph analysis in each function consistent, complete, and general for arbitrary networks. 
	May want to just have each function output a trajectory of population averages and sensitivities.
- Update the KMC analysis to be general for any reaction network (and read in the input file)
<<<<<<< HEAD
- Make separate functions for reading KMC output files and returning the information as variables
=======
- Post example output on google drive and share the link here.
	
Desired functionality
- Input reaction network in input files
- Choose between deterministic (ODE) and stochastic (KMC)
- Choose single time-scale (STS) or two time-scale mode (TTS)
- Specify which reactions are fast and slow
- Choose simulation parameters (time horizon, sampling times, number or replicates, etc.)
- Averageing over replicates for transient (don't do macro steady-state because that's more complicated)
- Simple post-processing scripts which read the output files and output trajectories and sensitivity estimates with error estimates
>>>>>>> origin/master

Future Work
- Port the STS and TTS KMC implementations into Fortran for speed