# signet: Selection Inference in Gene Networks

## Introduction

`signet` is an R package to detect natural selection in gene networks.

We will call gene network every type of data involving interactions
between genes or proteins. For example, it can be a gene regulatory network, 
a biological pathway, or a protein-protein interactions dataset.

Pritchard & Di Rienzo (2010) argued that many, or most, adaptive events in 
natural populations occur by the evolution of polygenic traits, rather than 
via the fixation of single beneficial mutations . Recent genome-wide association
studies in humans (Yang et al. 2010; Stranger et al. 2011) and in various model 
organisms (Atwell et al. 2010; Jumbo-Lucioni et al. 2010) have indeed  confirmed
that  variation  at  many  important traits is controlled by a large number of 
loci dispersed throughout the genome. Polygenic adaptation typically involves 
small allele frequency changes at many loci (Mackay et al. 2009). Therefore, 
quantitative trait loci (QTL) involved in polygenic adaptation will go 
undetected using methods for detecting the molecular signature of selective 
sweeps at individual loci (Pritchard & Di Rienzo 2010).

Here, we consider the information given by the network information as 
a prior concerning the target of polygenic selection. Indeed, it is more likely
that genes together under selection are found in the same biological network
because they are all related to the same phenotype.

The methodology implemented in the signet package is then an extension of 
genome scans for selection to gene networks. Using a statistic to measure
selection, the idea is to find in gene networks high scoring subgroups of genes.

## Workflow

### General principle

The general idea is simple: we have a gene list with prealably defined scores 
(e.g. a differentiation measure like the Fst) and we want to find gene networks
yielding a highest score than expected under a null hypothesis.

To do so, we will use an interaction database (e.g. biological pathways) and
search in these data for high-scoring subnetworks. Then, 

<img src="misc/workflow.png" width="400">

### Search algorithm

The method implemented is based on Ideker et al. (2002) heuristics, 
but several improvements are considered. You can see below an animation 
representing a run of the simulated annealing algorithm used in the package.

![simulatedAnnealing](misc/anim_50fps.gif)

As you can see, as we add or remove new genes in the active subnetwork (in red),
the score is maximized as we iterate.

You can read a more detailed description of the algorithm by clicking 
[here](misc/methodo.md).

### Testing the significance of the subnetworks scores

As the high-scoring search procedure tends to bias the p-values distribution 
towards low values, the test will be non-conservative if we use the background 
distribution to compute the significance of the subnetwork score. Therefore, 
the test implemented uses this high-scoring search procedure. For N iterations, 
gene scores are permuted. Then, a pathway is randomly sampled (the sampling 
probability being conditioned by the number of genes in the pathways) and 
the search algorithm is applied to this pathway. The score of the high-scoring 
subnetwork found in the randomized data is computed and its distribution is
generated. This is the null distribution of the test.

Finally, a correction for multiple testing is highly recommended as you usually
apply this test to hundreds of pathways. We advise to use the FDR 
method of X et al. (XXX) implemented in the R package `qvalue`.

## A walkthrough example

### Installation

There is no official release of the `signet` package at the moment. 
But you can install the development version on GitHub using the `devtools` 
package (`Rtools` must also be installed and properly configured):

```r
#install.packages('devtools')
devtools::install_github('algorythmes/signet')
```

### Data

We will use KEGG Pathways data, and genetic data from Daub et al. (2013), 
consisting in corrected Fst computed over 53 human populations, 
for more than 17,000 genes.

```r
data(keggPathways);data(zScores)
```

### Analysis

First, we generate the background distribution of the subnetworks scores 
for subnetworks of size f from 1 to 200 (the size of the biggest KEGG pathway). 
This may be a little long, so you can use `data(nullDistExample)` instead.

```r
nullDistribution(keggPathways,zScores,iterations = 10000)
```
Then, we apply the simulated annealing algorithm 
on pathways of your choice. Pathways must be in the `graphNEL` format. 
You can provide the `searchSubnet()` function a graph list, or a single graph.

```r
searchSubnet(keggPathways[[1]],zScores,iterations = 10000)
```

## References

Daub, J. T., Hofer, T., Cutivet, E., Dupanloup, I., Quintana-Murci, L., 
Robinson-Rechavi, M., & Excoffier, L. (2013). Evidence for polygenic 
adaptation to pathogens in the human genome. Molecular biology and evolution, 
30 (7): 1544-1558.

Sales, G., Calura, E., Cavalieri, D., & Romualdi, C. (2012). 
graphite-a Bioconductor package to convert pathway topology to gene network. 
BMC bioinformatics, 13(1), 20.
