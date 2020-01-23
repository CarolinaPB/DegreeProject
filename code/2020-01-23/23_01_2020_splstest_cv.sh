#!/bin/bash -l
#SBATCH -A g2019026
#SBATCH -p core #partition
#SBATCH -n 10 #number of tasks
#SBATCH -t 05:00:00 #total run time of the job allocation
#SBATCH --output="%A_%a.out"
#SBATCH --error="%A_%a.error"
#SBATCH --mail-type=ALL
#SBATCH --mail-user carolina.pitabarros@gmail.com

Rscript /home/carolpb/DegreeProject/code/2020-01-23/23_01_2020_spla_test.R