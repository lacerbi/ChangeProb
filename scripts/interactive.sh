#!/bin/bash
qsub  -lmem=4096mb -lnodes=1:ppn=1 -lwalltime=4:00:00 -qinteractive -I  -M${USER}@nyu.edu
