# signet: Selection Inference in Gene Networks

## Introduction

`signet` is an R package to detect natural selection in gene networks.

We call gene network every type of data involving interactions
between genes or proteins. For example, it can be a gene regulatory network, 
a biological pathway, or a protein-protein interactions dataset.

Natural selection, as an evolutionary force, shapes patterns of genetic
diversity in natural populations. 
Complex traits, polygenic selection. 
Here, we consider the information given by the
network information as a prior concerning the target of polygenic selection.

The methodology implemented in the signet package is an extension of 
genome scans for selection to gene networks. Using a statistic to measure
selection.

## Workflow

### General principle

<img src="misc/workflow.png" width="100">
![workflow](misc/workflow.png)

### Search algorithm

The method implemented is based on Ideker et al. (2002) heuristics, 
but several improvements are considered. Here is an animation representing a 
run of the simulated annealing algorithm used in the package:

![simulatedAnnealing](misc/anim_50fps.gif)

You can learn more about the methodology by clicking [here](misc/methodo.md).

## Application

### Installation

There is no official release of the `signet` package at the moment. 
But you can install the development version on GitHub using the `devtools` 
package (`Rtools` must also be installed and properly configured):

```
#install.packages('devtools')
devtools::install_github('algorythmes/signet')
```

### Data

We will use KEGG Pathways data, and genetic data from Daub et al. (2013), 
consisting in corrected Fst computed over 53 human populations, 
for more than 17,000 genes.

```
data(keggPathways);data(zScores)
```

### Analysis

First, we generate the background distribution of the subnetworks scores 
for subnetworks of size f from 1 to 200 (the size of the biggest KEGG pathway). 
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

