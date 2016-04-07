#!/bin/csh
#
# LSF batch script to run the test MPI code
#
#BSUB -P UYLE0010                       # Project 99999999
#BSUB -a poe                            # select poe
#BSUB -x                                # exclusive use of node (not_shared)
#BSUB -n 1                             # number of total (MPI) tasks
#BSUB -R "span[ptile=1]"               # run a max of 32 tasks per node
#BSUB -J runscript                           # job name
#BSUB -o runarwpost%J.out                     # output filename
#BSUB -e runarwpost%J.err                     # error filename
#BSUB -W 4:00                           # wallclock time
#BSUB -q regular                       # queue
#
mpirun.lsf ./script_figures.sh >script_figures.out
exit
