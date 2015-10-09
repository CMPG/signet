# signet: Selection Inference in Gene Networks

`signet` is a R package designed to handle the evolutionary analysis of gene networks.

From a biological pathways database, and a gene scores list, one could search for high-scoring subnetworks in biological pathways.

The method implemented is based on Ideker et al. (2002) heuristics, but several improvements are considered.

## Installation

There is no official release of the signet package at the moment. But you can install the development version using the `devtools` package (Rtools must also be installed and properly configured):

```
#install.packages('devtools')
devtools::install_github('algorythmes/signet')
```

## Methodology

The principle is as follow: 

Here is an animation to look at a run of the algorithm:

![simulatedAnnealing](misc/anim_50fps.gif)

## Application

Data: we will use KEGG Pathways data, and scores used in Daub et al. (2013), consisting in scaled differentiation scores over 53 human populaions, for more than 17,000 genes.

```
data(keggPathways);data(zScores)
```

First, we generate the null distribution of the subnetworks scores for subnetworks of size kmin to kmax. This may be a little long, so you can use `data(nullDistExample)` here.

```
nullDistribution(keggPathways,zScores,iterations = 10000)
```
Then, you can apply the simulated annealing algorithm on pathways of your choice. You can provide the `searchSubnet()` function a graph list, or a single graph.

```
searchSubnet(keggPathways[[1]],zScores,iterations = 10000)
```

