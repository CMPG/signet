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
<<<<<<< HEAD

data(scores);data(kegg)

#Generate background distribution:
bkgd_dist <- backgroundDist(kegg_clean,scores)

#Run simulated annealing on the first 10 pathways
HSS <- searchSubnet(kegg_clean[1:10],scores,bkgd_dist)

#Generate the null distribution

#Results

=======
>>>>>>> 2cd560f40a476c0551af96785889d39bca6429f1
```
