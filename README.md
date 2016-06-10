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