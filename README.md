# signet: Selection Inference in Gene Networks

## Introduction

`signet` is a R package to detect natural selection in gene networks.

### What is a gene network?

A gene network. It can be a gene regulatory network, a biological pathway,
protein-protein interactions, or every kind of data involving interactions
between genes.

### What is selection in gene networks?

Natural selection is an evolutionary force shaping patterns of genetic
diversity. Polygenic selection. Here, we consider the information given by the
network information as a prior concerning the target of polygenic selection.

### How to detect selection in gene networks?

The methodology implemented in the signet package is an extension of 
genome scans for selection to gene networks. Using a statistic to measure
selection.

## Installation

There is no official release of the signet package at the moment. 
But you can install the development version using the `devtools` package 
(Rtools must also be installed and properly configured):

```
#install.packages('devtools')
devtools::install_github('algorythmes/signet')
```

## Methodology

The method implemented is based on Ideker et al. (2002) heuristic, 
but several improvements are considered. Here is an animation representing a 
run of the simulated annealing algorithm used in the package:

![simulatedAnnealing](misc/anim_50fps.gif)

You can learn more about the methodology by clicking [here](misc/methodo.md).

## Application

Data: we will use KEGG Pathways data, and genetic data from Daub et al. (2013), 
consisting in corrected Fst computed over 53 human populaions, 
for more than 17,000 genes.

```
data(keggPathways);data(zScores)
```

First, we generate the null distribution of the subnetworks scores 
for subnetworks of size kmin to kmax. 
This may be a little long, so you can use `data(nullDistExample)` instead.

```
nullDistribution(keggPathways,zScores,iterations = 10000)
```
Then, we apply the simulated annealing algorithm 
on pathways of your choice. Pathways must be in the `graphNEL` format. 
You can provide the `searchSubnet()` function a graph list, or a single graph.

```
searchSubnet(keggPathways[[1]],zScores,iterations = 10000)
```

