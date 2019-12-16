library("Hmisc") # for rcorr

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

cor_traits <- rcorr(as.matrix(phenotype, rownames = TRUE))

cor_traits.cor <- cor_traits$r
cor_traits.p <- cor_traits$P

cor_traits.cor[upper.tri(cor_traits.cor)] <- NA
cor_traits.p[upper.tri(cor_traits.p)] <- NA

my_cor_matrix <- flat_cor_mat(cor_traits$r, cor_traits$P)
head(my_cor_matrix)
colnames(my_cor_matrix) <- c("Trait1", "Trait2", "cor", "pval")
my_cor_matrix <- na.omit(my_cor_matrix)
head(my_cor_matrix)

my_cor_matrix <- data.table(my_cor_matrix)
fwrite(my_cor_matrix, "results/2019-12-13/genes_correlation.gz")
