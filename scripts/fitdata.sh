#!/bin/sh

module purge
#. /etc/profile.d/modules.sh

# Use Intel compiler
module load matlab/2016b
export MATLABPATH=${MATLABPATH}:${HOME}/${PROJECT}:${HOME}/MATLAB
source /home/${USER}/MATLAB/setpath.sh

#Check if running as an array job
if [[ ! -z "$PBS_ARRAYID" ]]; then
	IID=${PBS_ARRAYID}
fi
#Check if running as an array job
if [[ ! -z "$SGE_TASK_ID" ]]; then
        IID=${SGE_TASK_ID}
fi
if [[ ! -z "$SLURM_ARRAY_TASK_ID" ]]; then
        IID=${SLURM_ARRAY_TASK_ID}
fi

# Run the program
#echo ${WORKDIR} ${PROJECT} ${IID}.job

cat<<EOF | matlab -nodisplay
addpath(genpath('${HOME}/MATLAB'));
addpath(genpath('${HOME}/${PROJECT}'));
cd('${WORKDIR}');
changeprob_runfit(${IID},${FIXNOISE},${GRIDSIZE});
EOF
