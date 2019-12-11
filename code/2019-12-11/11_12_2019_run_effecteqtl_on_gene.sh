#!/bin/bash -l
#SBATCH -A g2019026
#SBATCH -p core #partition
#SBATCH -n 4 #number of tasks
#SBATCH -t 12:00:00 #total run time of the job allocation
#SBATCH --array=1-4
#SBATCH --output="%A_%a.out"
#SBATCH --error="%A_%a.error"
#SBATCH --mail-type=ALL
#SBATCH --mail-user carolina.pitabarros@gmail.com

Rscript /home/carolpb/DegreeProject/code/2019-12-09/09_12_2019_effecteqtl_on_genes.R $SLURM_ARRAY_TASK_ID
