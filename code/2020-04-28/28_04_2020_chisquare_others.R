## Chi-square test
H0: the categories are independent. 
H1: the categories are not independent

```{r}
table.chi2 <- table(fract_less10_toplot$fraction, fract_less10_toplot$causal)
chi2.result <- chisq.test(table.chi2)
print(chi2.result)

if (chi2.result$p.value < 0.01){
  print(paste0("since the p-value (",round(chi2.result$p.value, digits = 4), ") is less than the significance level 0.01, we reject the null hypothesis - there's an association between the causal categories")) 
} else {
  print("We don't reject the null hypothesis: there is no association between the causal categories")
}
```

Testing the difference between "causal in hotspot" and "not causal not in hotspot"
```{r}
table.chi2 <- table(fract_less10_toplot[causal %in% c("causal in hotspot", "not causal not in hotspot")]$fraction, fract_less10_toplot[causal %in% c("causal in hotspot", "not causal not in hotspot")]$causal)
chi2.result <- chisq.test(table.chi2)
print(chi2.result)

if (chi2.result$p.value < 0.01){
  print(paste0("since the p-value (",round(chi2.result$p.value, digits = 4), ") is less than the significance level 0.01, we reject the null hypothesis - there's an association between the causal categories")) 
} else {
  print("We don't reject the null hypothesis: there is no association between the causal categories")
}
```