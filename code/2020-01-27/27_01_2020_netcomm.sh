#!/bin/bash -l
#SBATCH -A snic2019-8-367
#SBATCH -p core #partition
#SBATCH -n 3 #number of tasks
#SBATCH -t 01:00:00 #total run time of the job allocation
#SBATCH --output="%A_%a.out"
#SBATCH --error="%A_%a.error"
#SBATCH --mail-type=ALL
#SBATCH --mail-user carolina.pitabarros@gmail.com

Rscript /home/carolpb/DegreeProject/code/2020-01-27/27_01_2020_netcomm.R
