#!/bin/bash
PROJECT="ChangeProb"
SHORTNAME="CP"
BASEDIR="/home/${USER}/${PROJECT}"
SOURCEDIR="${BASEDIR}"
JOBSCRIPT="${BASEDIR}/scripts/fitdata.sh"

#Job parameters
RUN=${1}
WORKDIR="/scratch/${USER}/${PROJECT}/run${RUN}"
mkdir ${WORKDIR}
cd ${WORKDIR}
MAXID=290
RUNTIME=48:00:00
FIXNOISE="[]"
GRIDSIZE="50"

#Job list is second argument
if [[ ! -z "$2" ]]; then
        JOBLIST=$2
else
	#Get all files in directory
	JOBLIST="1-${MAXID}"
fi

#RESOURCES="nodes=1:ppn=1,mem=6GB,walltime=${RUNTIME}"
RESOURCES="nodes=1:ppn=1,mem=3GB,walltime=${RUNTIME}"

#Convert from spaces to commas
JOBLIST=${JOBLIST// /,}
echo JOBS $JOBLIST

JOBNAME=${SHORTNAME}${RUN}
qsub -t ${JOBLIST} -v PROJECT=${PROJECT},RUN=${RUN},MAXID=$MAXID,WORKDIR=$WORKDIR,USER=$USER,FIXNOISE=${FIXNOISE},GRIDSIZE=${GRIDSIZE} -l ${RESOURCES} -M ${USER}@nyu.edu -N ${JOBNAME} ${JOBSCRIPT}
