#!/bin/sh
#PBS -o localhost:${PBS_O_WORKDIR}/
#PBS -e localhost:${PBS_O_WORKDIR}/
#PBS -M ehn222@nyu.edu
#PBS -q normal

module purge
#. /etc/profile.d/modules.sh

# Use Intel compiler
module load matlab
export MATLABPATH=${MATLABPATH}:/home/${USER}/${PROJECT}:/home/${USER}/MATLAB
source /home/${USER}/MATLAB/setpath.sh

#Check if running as an array job
if [[ ! -z "$PBS_ARRAYID" ]]; then
	IID=${PBS_ARRAYID}
fi
#Check if running as an array job
if [[ ! -z "$SGE_TASK_ID" ]]; then
        IID=${SGE_TASK_ID}
fi

# Run the program
#echo ${WORKDIR} ${PROJECT} ${IID}.job

cat<<EOF | matlab -nodisplay
addpath(genpath('/home/${USER}/MATLAB'));
addpath(genpath('/home/${USER}/${PROJECT}'));
cd('${WORKDIR}');
changeprob_runfit(${IID});
EOF
