#PBS -S /bin/bash
#PBS -l nodes=10,tpn=10,naccesspolicy=shared
#PBS -l walltime=00:20:00
#PBS -N log_NEB
#PBS -t 0-3
#PBS -o logs/log_NEB_forces

##########################
dirs=(
im_0.00
im_1.00
)
#dirs=( $(find -maxdepth 1 -type d  -name "im_*" -printf '%P\n' ) )

####PBS -t 1-${#dirs[@]}
##########################

# load modules
module purge
module load null torque moab
module load intel/parallelxe/2016

ulimit -s unlimited
ulimit -c unlimited
export OMP_STACKSIZE=64M

# Construct directory name for current PBS array job index
dir=${dirs[$PBS_ARRAYID]}

# Change to directory from which we want to run the input file
cd $PBS_O_WORKDIR/$dir

echo FILES:
ls

# Set up threading and machinefile
cat $PBS_NODEFILE | uniq > nodes_$PBS_ARRAYID
ncores=`wc -l $PBS_NODEFILE | awk '{ print $1 }'`
nthreads=`grep -m1 threads_max $PBS_O_WORKDIR/$dir/*.dat |  awk '{ print $3 }'`
nprocs=`echo $(($ncores / $nthreads))`

export OMP_NUM_THREADS=$nthreads

# Output execution information to output file
touch onetep.out
cat nodes_$PBS_ARRAYID >> onetep.out
echo "Processes: " $nprocs "Threads: " $nthreads "Cores: " $ncores >> onetep.out
echo "mpirun -n $nprocs -machinefile nodes_$PBS_ARRAYID onetep.iridis4.intel16.omp.scalapack onetep.dat 2> onetep.err" >> onetep.out

# Run calculation
mpirun -np $nprocs -machinefile nodes_$PBS_ARRAYID -bootstrap rsh ../hpc_files/onetep.iridis4.intel16.omp.scalapack *.dat >> onetep.out 2> onetep.err



