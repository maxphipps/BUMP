#/bin/bash

# PBS submission script submitter
# -> Submits a series of NEB optimisation steps.
#
# Max Phipps, 2017

#NUM_ITER=3


# Ensure log/ directory exists
mkdir -p logs

# Run the array of intermediate replica singlepoint calculations
# to calculate the current forces 
JOB_ID_1=`qsub sub_neb_forces.sh`

# Once all replica jobs are complete, 
# iterate the nudged elastic band optimisation by steepest descent
# (using the above calculated forces)
JOB_ID_2=`qsub -W depend=afterokarray:$JOB_ID_1 sub_neb_iteration.sh`

for i in {1..3}
#START=1
#for (( i=$START; i<=$END; i++ ))
do
   echo "NEB iteration: $i"

   # Once the NEB steepest descent iteration python script is complete,
   # Run an array of intermediate replica singlepoint calculations
   # to calculate the forces 
   JOB_ID_1=`qsub -W depend=afterok:$JOB_ID_2 sub_neb_forces.sh`
   
   # Once all replica jobs are complete, 
   # iterate the nudged elastic band optimisation by steepest descent
   # (using the above calculated forces)
   JOB_ID_2=`qsub -W depend=afterokarray:$JOB_ID_1 sub_neb_iteration.sh`

done
