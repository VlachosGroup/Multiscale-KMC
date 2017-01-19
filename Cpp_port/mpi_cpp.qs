# Squidward C++ MPI parallelization submit file

#!/bin/bash
#$ -cwd
#$ -j y
#$ -N parallel_cpp
#$ -S /bin/bash
#$ -pe openmpi-smp 16
#

# Get our environment setup:

# The  executable:
mpiexec -n 5 ./a.out 