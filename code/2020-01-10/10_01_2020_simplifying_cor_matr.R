library(data.table)

path <- "/home/carolpb/DegreeProject/"

flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
  # Column 1 : row names (variable 1 for the correlation test)
  # Column 2 : column names (variable 2 for the correlation test)
  # Column 3 : the correlation coefficients
  # Column 4 : the p-values of the correlations
  
  library("dplyr")
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
}

load(file = paste0(path, "results/2020-01-08/full_cor.Rdata"))

cor_traits.cor <- cor_traits$r
cor_traits.p <- cor_traits$P

my_cor_matr_flat <- flat_cor_mat(cor_traits.cor,cor_traits.p)
my_cor_matr_flat <- data.table(my_cor_matr_flat)
cor_matr <- my_cor_matr_flat[!duplicated(t(apply(my_cor_matr_flat, 1, sort))), ]

fwrite(cor_matr, paste0(path, "results/2020-01-10/cor_matr_unique.gz"))