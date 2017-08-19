#PBS -S /bin/bash
#PBS -l nodes=1:ppn=1,naccesspolicy=shared,pmem=50mb
#PBS -l walltime=00:00:20
#PBS -N log_NEB_iter
#PBS -o logs/log_NEB_iter
#PBS -q test

##PBS -l nodes=1,tpn=1,naccesspolicy=shared
##PBS -l nodes=1:ppn=1,naccesspolicy=shared

# load modules
module purge
module load null torque moab
module load python/2.7.5

# Change to directory submission script was submitted from
cd $PBS_O_WORKDIR

# Iterate the NEB optimisation
python2.7 neb.py -iter

