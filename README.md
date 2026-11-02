# sLDtest

The 'sLDtest' R package implements a hypothesis testing method of structural linkage disequilibrium (sLD), as well as a forward-generating probability model for genotypes in sLD.

## Installation 

```R
install.packages("pak")
pak::pak("StoreyLab/sLDtest")
```

## Examples 

### Input genotype data 

'sLDtest' is designed to be applied to genotype data like in the example below. These genotypes were simulated under an admixture model with $K=4$ admixed populations using the package 'bnpsd'.

```R
# Load input genotype data matrix (I loci by J individuals)
data(genodata)
```

### Run the sLD test 

The implementation of the sLD test requires only genotypes to be supplied by the user. However, individual-specific allele frequencies (IAFs) must be estimated from these genotypes using any one of a wide range of existing methods. In our analyses, we use the package 'alstructure' to estimate IAFs. 

```R
# Estimate IAFs 
X <- genodata
J_ind <- ncol(X) # get the number of individuals 
Pi_hat <- alstructure::estimate_F(X, d = 4)$F_hat # assume K known for simplicity
Pi_hat[Pi_hat < 1/(2*J_ind)] <- 1/(2*J_ind) # truncate estimates
Pi_hat[Pi_hat > 1 - 1/(2*J_ind)] <- 1 - 1/(2*J_ind)
```

After estimating IAFs, we can compute sLD test statistics. The user supplies a list of locus pairs to perform the sLD test on as follows. Here, we simply randomly select 1000 pairs of loci. 

```R
n_pairs <- 1000
pair_ids <- t(replicate(n_pairs, sample(nrow(X), size = 2, replace = FALSE)))
test_stats <- sLD(X, Pi_hat, pair_ids)
``` 

We then compute $p$-values as follows. We use a Normal(0,1) null distribution by default, though we can also supply an empirical null distribution if necessary.

```R
p_vals <- compute_p_val(test_stats)
```

To reproduce the results in the paper "A test for linkage disequilibrium in structured populations" by Alvin J. Zhang and John D. Storey, please visit the following repository:  
  
https://github.com/StoreyLab/sLD-manuscript-analysis

## Citations 

Alvin J. Zhang, John D. Storey. 2026. "A test for linkage disequilibrium in structured populations." bioRxiv [doi:10.64898/2026.10.18.759967](https://doi.org/10.64898/2026.10.18.759967) 2026-10-19. 
For users of 'alstructure', please cite: Irineo Cabreros, John D. Storey. 2019. "A likelihood-free estimator of population structure bridging admixture models and principal components analysis." doi: [10.1534/genetics.119.302159](https://doi.org/10.1534/genetics.119.302159) 2019-08-01. 
For users of 'bnpsd', please cite: Alejandro Ochoa, John D Storey.  2021.  "Estimating FST and kinship for arbitrary population structures." PLoS Genet 17(1): e1009241. PubMed ID 33465078. [doi:10.1371/journal.pgen.1009241](https://doi.org/10.1371/journal.pgen.1009241). bioRxiv [doi:10.1101/083923](https://doi.org/10.1101/083923) 2016-10-27.

 