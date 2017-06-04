# signet: Selection Inference in Gene NETworks

We have developed a method to detect selection in biological pathways. 
The general idea is to search for subsets of interacting genes within biological pathways
that present unusual features. This search 
is a typical combinatorial optimization problem that can be solved using a 
heuristic approach like simulated annealing. We implemented such an algorithm 
to search for high-scoring subnetworks of genes in biological pathways data. 
We also developed a significance test procedure that explicitly takes into 
account the optimization process involved in this search.

This workflow is implemented in the `signet` package presented here. It is an 
extension of genome scans for selection to gene networks. Using a statistic to 
measure selection, the idea is to find high scoring subnetworks 
of genes within biological pathways.

## General principle

The general idea is simple: we consider a gene list with prealably defined scores 
(e.g. a differentiation measure like the Fst) and we want to find gene networks
presenting a score higher than expected under the null hypothesis.

To do so, we will use biological pathways databases converted as gene networks 
and search in these graphs for high-scoring subnetworks.

Details about the algorithm can be found <a href="http://biorxiv.org/content/early/2017/04/18/128306">here</a>.

The method implemented is based on Ideker et al. (2002) heuristics, 
but several improvements are considered. You can see below an animation 
representing a run of the simulated annealing algorithm used in the package.

<p align="center"><img src="misc/anim_50fps.gif"></p>
<p align="center">Simulated annealing run in a gene network.</p>

## Walkthrough example

### Installation

There is no official release of the `signet` package at the moment. 
But you can install the development version on GitHub using the `devtools` 
package (`Rtools` must also be installed and properly configured):

```{r}
#install.packages('devtools')
devtools::install_github('CMPG/signet')
library(signet)
```
<<<<<<< HEAD
=======
<<<<<<< HEAD

### Input

We will use gene scores obtained from Daub et al. (2013) and run the search
in KEGG pathways.

```{r}
library(graphite)

kegg <- pathways("hsapiens", "kegg") # KEGG pathways
kegg <- lapply(kegg,pathwaysGraph) # pathways objects need to be converted as graphs
data(daub2013) # gene scores
```

### Analysis

NB: This procedure is still a bit slow (but doable in a few hours on
a single computer); a C++ implementation will be soon considered.

First, we need to generate the background distribution of the subnetworks scores 
for all possible subnetwork sizes.
You can skip this step and run `data(backgroundDist)` instead.

```{r}
bkgd <- backgroundDist(kegg,scores)
```

Then, we apply the simulated annealing algorithm on a list of pathways
(here, the first 10 KEGG pathways). Pathways must be in the `graphNEL` format.

```{r}
analysis <- searchSubnet(pathways = kegg[1:10],
                             score = scores,
                             background = bkgd)
```

This function returns a list of N signet objects (corresponding to N pathways),
including a table with the whole list of genes, their scores, 
and a boolean indicating if they are found in the high-scoring subnetwork.

To assess the significance of the subnetwork scores, we need to generate an 
empirical null distribution.

```{r}
null <- nullDist(kegg, 1000)
# and compute p-values:
### 
```

### Extracting the results

You can generate a table using the `summary()` function and then write the 
results in your working directory.

```{r}
tab <- summary(analysis)
write.table(tab, sep="\t", file="results.tsv")
```

Plot the results in R

To get a better representation of the networks, we advise to use the Cytoscape software.
The package RCytoscape
Works only with Cytoscape v.X.X.X

Plotting the results with Cytoscape.
Single pathway.

Merged significant pathways.
=======
>>>>>>> b488d6d1fb3c0e223e815a75e3374621decb452f
>>>>>>> 6901074f9e2d9e453dc937b9474634192e77ab58
