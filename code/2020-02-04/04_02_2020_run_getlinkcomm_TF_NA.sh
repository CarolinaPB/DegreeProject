#!/bin/bash -l
#SBATCH -A g2019026
#SBATCH -p core #partition
#SBATCH -n 4 #number of tasks
#SBATCH -t 24:00:00 #total run time of the job allocation
#SBATCH --output="%A_%a.out"
#SBATCH --error="%A_%a.error"
#SBATCH --mail-type=ALL
#SBATCH --mail-user carolina.pitabarros@gmail.com

Rscript /home/carolpb/DegreeProject/code/2020-02-04/04_02_2020_getlinkcomm_TF_TNA.R