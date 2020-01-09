#!/bin/bash -l
#SBATCH -A g2019026
#SBATCH -p core #partition
#SBATCH -n 8 #number of tasks
#SBATCH --output="%A_%a.out"
#SBATCH --error="%A_%a.error"
#SBATCH --qos=short # for test jobs
#SBATCH --mail-type=ALL
#SBATCH --mail-user carolina.pitabarros@gmail.com

Rscript /home/carolpb/DegreeProject/code/2020-01-09/09_01_2020_anova_uppmax_parallel.R