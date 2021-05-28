#!/bin/sh
########################################
#                                      #
# GE job script for ECDF Cluster       #
#                                      #
# Template by ECDF System Team         #
# ecdf-systems-team@lists.ed.ac.uk     #
# Specific modifications made by       #
# A. A. Bloom (UoE, now at JPL)        #
#                                      #
########################################

# Grid Engine options

#$ -cwd
###$ -P geos_gcri_cardamom

# Initialise environment module
. /etc/profile.d/modules.sh

# Use Intel compiler

module load intel

#THIS SCRIPT MUST BE ACCOMPANIED BY CARDAMOM_ECDF_EXECUTABLES_LIST.txt IN THE SAME DIRECTORY
# arguments are start and end lines!

#for line_number in $(seq $2 1 $3)
#do
#task=$( cat $1CARDAMOM_ECDF_EXECUTABLES_LIST.txt | sed ${line_number}\!d )
task=$( cat $1CARDAMOM_ECDF_EXECUTABLES_LIST.txt | sed $SGE_TASK_ID\!d )
command ${task}

#done
