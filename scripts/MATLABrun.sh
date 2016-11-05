#!/bin/bash
module purge
module load matlab
export MATLABPATH=${MATLABPATH}:/home/${USER}/MATLAB
matlab -nodisplay
