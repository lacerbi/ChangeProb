#!/bin/bash
PROJECT="ChangeProb"
SHORTNAME="CP"
BASEDIR="${HOME}/${PROJECT}"
SOURCEDIR="${BASEDIR}"
JOBSCRIPT="${BASEDIR}/scripts/fitdata.sh"

#Job parameters
RUN=${1}
WORKDIR="${SCRATCH}/${PROJECT}/run${RUN}"
mkdir ${WORKDIR}
cd ${WORKDIR}
MAXID=290
RUNTIME=24:00:00
FIXNOISE="[]"
GRIDSIZE="[]"

#Job list is second argument
if [[ ! -z "$2" ]]; then
        JOBLIST=$2
else
	#Get all files in directory
	JOBLIST="1-${MAXID}"
fi

NODES="1"
PPN="1"
MEM="3GB"
RESOURCES="nodes=${NODES}:ppn=${PPN},mem=${MEM},walltime=${RUNTIME}"

#Convert from spaces to commas
JOBLIST=${JOBLIST// /,}
echo JOBS $JOBLIST

JOBNAME=${SHORTNAME}${RUN}

if [ ${CLUSTER} = "Prince" ]; then
        # running on Prince
        sbatch --verbose --array=${JOBLIST} --mail-type=FAIL --mail-user=${USER}@nyu.edu --mem=${MEM} --time=${RUNTIME} --nodes=${NODES} --ntasks-per-node=${PPN} --export=PROJECT=${PROJECT},RUN=${RUN},MAXID=$MAXID,WORKDIR=$WORKDIR,USER=$USER,FIXNOISE=${FIXNOISE},GRIDSIZE=${GRIDSIZE} --job-name=${JOBNAME} ${JOBSCRIPT}
else
	qsub -t ${JOBLIST} -q normal -v PROJECT=${PROJECT},RUN=${RUN},MAXID=$MAXID,WORKDIR=$WORKDIR,USER=$USER,FIXNOISE=${FIXNOISE},GRIDSIZE=${GRIDSIZE} -l ${RESOURCES} -M ${USER}@nyu.edu -N ${JOBNAME} ${JOBSCRIPT}
fi
