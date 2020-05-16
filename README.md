# signet: Selection Inference in Gene NETworks

The `signet` R package implements a method to detect selection in biological 
pathways. The general idea is to search for gene subnetworks within biological 
pathways that present unusual features, using a heuristic approach 
(simulated annealing). Details about the algorithm can be found in
<a href="https://doi.org/10.1093/nar/gkx626">Gouy et al. (2017)</a>.

### Installation

`signet` is available on Bioconductor.

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("signet")
```

You can also install the development version:

```{r}
## install devtools and signet dependencies
# install.packages('devtools')
# if (!requireNamespace("BiocManager", quietly=TRUE))
    # install.packages("BiocManager")
# BiocManager::install("graph")
# BiocManager::install("RBGL")
# install.packages("igraph")

devtools::install_github('CMPG/signet')
library(signet)
```

### Usage

You can browse `signet` vignette to learn how to use the package with 
a walkthrough example.

```{r}
library(signet)
vignette("signet")
```

### Citation

Please cite <a href="https://doi.org/10.1093/nar/gkx626">this paper</a> if you use `signet` for your project:

* Gouy A, Daub JT & Excoffier L (2017) Detecting gene subnetworks under 
selection in biological pathways. Nucleic Acids Research 45(16): e149
