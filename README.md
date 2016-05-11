# Multiscale-KMC-with-SA
Two time-scale well-mixed kinetic Monte Carlo. Multiscale likelihood ratio sensitivity analysis is done on the fly. Implemented in Matlab. Single time-scale and ODE analogues are also included.

**Citations**  
See license or references page in the wiki.

**Workflow**  
Please see wiki page for detailed instructions.
- Input reaction network in input files
- Specify which reactions are fast and slow
- Choose numerical parameters (time horizon, sampling times, number or replicates, etc.)
- Choose a model: deterministic (ODE) vs. stochastic (KMC), single time-scale (STS) vs. two time-scale mode (TTS)
- Averageing over replicates for transient (don't do macro steady-state because that's more complicated)
- Simple post-processing scripts which read the output files and output trajectories and sensitivity estimates with error estimates

**To Do**  
- Read reaction network and rate information from an input file.
- Implement analysis for KMC STS

- Error handle: if there is no equilibrium with choice of microscale reactions (or there are none)

- Include a latex document with a simple example worked out (with data). There is an example in Research_Files/Thesis/MultiscaleSA
<<<<<<< HEAD
- Post example output on google drive and share the link here.
- Create sample data for KMC STS
- Post example output on google drive and share the link here.

Future Work
=======
- Update the KMC analysis to be general for any reaction network (and read in the input file)
- Post example output on google drive and share the link here.
- Create sample data for KMC STS
>>>>>>> origin/master
- Port the STS and TTS KMC implementations into Fortran for speed

