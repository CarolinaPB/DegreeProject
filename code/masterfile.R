# Project code Johan
# g2019026

# To use Johan's project code, the batch jobs need to be sent to snowy https://www.uppmax.uu.se/support/user-guides/snowy-user-guide/
# SNOWY
squeue -M snowy
sbatch -M snowy slurm_script_file

# to check usage while job is running --> snowy
less /sw/share/slurm/snowy/uppmax_jobstats/*/<job_ID>
# to check usage while job is running --> snowy
less /sw/share/slurm/rackham/uppmax_jobstats/*/<job_ID>
  
# Other job stats
scontrol show job <job id>
scontrol -M snowy show job <job id>
  
  
#   ____________________________________________________________________________
#   Starting data                                                           ####
library(data.table)
phenotype <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/data/SI_Data_01_expressionValues.txt")
genotype <- fread("data/SI_Data_03_genotypes.txt")
eqtl_results <- fread("data/SI_Data_04_eQTL.csv")


#   ____________________________________________________________________________
#   General info                                                            ####

# expression levels in units of log2(TPM) for all genes and segregants

# 1012 haploid segregants (rows) X 42,052 markers (columns)
# -1 indicates BY allele, +1 indicates RM allele
# BY (i.e. reference) alleles are denoted by ‘-1’. RM alleles are denoted by ‘1’
# column names indicate chromosome:position_BYvariant/RMvariant
# postitions are based on the S.Cerevisiae SacCer3 genome build

# LOD stands for "logarithm of the odds." 
# In genetics, the LOD score is a statistical estimate of whether two genes, or a gene and a disease gene,
# are likely to be located near each other on a chromosome and are therefore likely to be inherited. 
# A LOD score of 3 or higher is generally understood to mean that two genes are located close to each other on the chromosome. 
# In terms of significance, a LOD score of 3 means the odds are a thousand to one that the two genes are linked, and therefore inherited together.

  
  
