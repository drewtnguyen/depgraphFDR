This is the `R` package that implements the procedures 
in INSERT URL HERE.

Install in `R` 
```{r}
devtools::install_github("drewtnguyen/depgraphFDR")
```

This package is under active development. Its API is subject to change. 

# Basic usage

We try our method with the running example graph from the paper. 

```{r}
library(depgraphFDR)
adjlist <- list(
  "1" = c(1, 2, 3),  
  "2" = c(1, 2, 3),  
  "3" = c(1, 2, 3, 4, 5),  
  "4" = c(3, 4),  
  "5" = c(3, 5)   
)

pvals = c(0.02, 0.02, 0.01, 0.02, 0.04)
rej_ind = IndBH(alpha = 0.05, pvals = pvals, adjlist = adjlist)
rej_ind2 = IndBH_plus(alpha = 0.05, pvals = pvals, adjlist = adjlist, recurse = 1)
```

