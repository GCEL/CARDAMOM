#!/bin/sh
########################################
#                                      #
# GE job script for ECDF Cluster       #
#                                      #
# Template by ECDF System Team         #
# ecdf-systems-team@lists.ed.ac.uk     #
# Specific modifications made by       #
# A. A. Bloom (UoE, now at JPL)        #
# T. L. Smallman (UoE)                 #
#                                      #
########################################

# Grid Engine options

#$ -cwd
###$ -P geos_gcri_cardamom

# Initialise environment module
. /etc/profile.d/modules.sh

# Use Intel compiler

module load intel
module load gcc

# THIS SCRIPT MUST BE ACCOMPANIED BY CARDAMOM_ECDF_EXECUTABLES_LIST.txt IN THE SAME DIRECTORY
# arguments are start and end lines!

task=$( cat $1CARDAMOM_ECDF_EXECUTABLES_LIST.txt | sed $SGE_TASK_ID\!d )
export OMP_NUM_THREADS=$NSLOTS
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_DYNAMIC=false
export OMP_STACKSIZE=1G
ulimit -s unlimited
command ${task}

