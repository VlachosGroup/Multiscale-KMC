#!/bin/bash
#$ -cwd
#$ -N MSA-KMC
#$ -pe openmpi-smp 16
#$ -q *@@3rd_gen
#$ -j y
#$ -S /bin/bash
#$ -o MSA-KMC.out

# squidward.che.udel.edu submission script for C++ with MPI parallelization

# Setup environment:
source /etc/profile.d/valet.sh
vpkg_require openmpi/1.6.3-gcc

# Begin the run by printing some basic info and then
# invoke mpiexec:
echo "GridEngine parameters:"	
echo "  nhosts         = $NHOSTS" 
echo "  nproc          = $NSLOTS" 
echo "  mpiexec        =" `which mpiexec` 
echo "  pe_hostfile    = $PE_HOSTFILE" 	
echo 
cat $PE_HOSTFILE 
echo 
echo "-- begin OPENMPI run --"
time mpiexec --n $NSLOTS ./MSA-KMC-MPI.x
echo "-- end OPENMPI run --"