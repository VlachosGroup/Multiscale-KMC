#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe mpich 1
#$ -N MSA-KMC
#$ -o MSA-KMC.out

# squidward.che.udel.edu submission script for a serial job
source /etc/profile.d/valet.sh
vpkg_require openmpi/1.6.3-gcc

# The  executable:
export KMC_EXE="./MSA-KMC.x"

# Simple summary:
echo ""
echo "Running on ${HOSTNAME} with job id ${JOB_ID}"
echo ""

time ${KMC_EXE}