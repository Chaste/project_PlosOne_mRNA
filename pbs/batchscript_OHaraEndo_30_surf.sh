#!/bin/bash

#PBS -V
#PBS -l select=1
#PBS -l walltime=01:00:00
#PBS -N Popn_30surf

# /home/blanca-rodriguez/walmsley/bin/toRunChaste

cd $PBS_O_WORKDIR

TAG=_30_surf

EXECUTABLE_PATH=/home/blanca-rodriguez/walmsley/soft/install_chaste/Chaste/projects/JohnW/build/intel/Sensitivity_Analysis
CHASTE_EXECUTABLE=TestSensitivityAnalysisOHaraEndoRunner
INPUT=/home/blanca-rodriguez/walmsley/soft/install_chaste/Chaste/projects/JohnW/test/data/Sensitivity/ChasteInput/input${TAG}/exp_design${TAG}
# write to time file
date > model.time

# launch processes
for INC in {0..7}; do
  let "id=$[PBS_ARRAY_INDEX]+$[INC]"
  $EXECUTABLE_PATH/${CHASTE_EXECUTABLE} --file ${INPUT}_${id}.dat --tag $TAG &
done

# wait for processes test1 to finish
sleep 4
while [ "$(pidof $CHASTE_EXECUTABLE)" ]; do
    sleep 60
done

# append to time file
date >> model.time
