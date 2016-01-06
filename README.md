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

### Input preparation

#### Biological networks conversion to graphs

In order to get information about relationships between genes,
we advise to use the Bioconductor package `graphite`, which implements a 
procedure to convert biological pathways data from different databases to 
graphs (for more details, see Sales et al., 2012). 
These databases have emerged as references in systems biology:
KEGG, Reactome, BioCarta, NCI, Panther and HumanCyc.

#### Correction for overlapping between pathways

We will first exclude all pathways completely included in at least another 
pathway (overlap = 1). I still need to implement a procedure to merge pathways 
with an overlapping level greater than a given threshold.

Finally, pathways with no connected component of size greater than 10 
(arbitratry threshold) are removed.

### Search algorithm

The method implemented is based on Ideker et al. (2002) heuristics, 
but several improvements are considered. You can see below an animation 
representing a run of the simulated annealing algorithm used in the package.

![simulatedAnnealing](misc/anim_50fps.gif)

As you can see, as we add or remove new genes in the active subnetwork (in red),
the score is maximized as we iterate.

#### Background distribution

To apply the algorithm to a list of pathways, you have to first generate the 
background distribution of your summary statistic.

The input is the list of all the pathways you want to analyze, and a list of 
scores for each gene. The background distribution or the statistic will be 
generated for each possible subnetwork size (kmin to kmax).

A pathway is randomly sampled (the sampling probability being conditioned by 
the number of genes in the pathway) and a connected subnetwork of size k is 
randomly picked (one gene is sampled in the pathway, then k genes in the 
boundary). The score of the subnetwork is then computed. This is done N times 
for each k, to get the background distribution of subnetworks scores, we just 
keep the mean and standard error of this distribution.

#### Algorithm for high-scoring subnetwork search

We consider that genes can yield two states: active or inactive.

The selected gene is not an articulation point of the subgraph, i.e. its 
removal doesn’t disconnect the active subgraph.

The selected gene is randomly picked in the following nodes: i) nodes in 
the boundary; ii) leaves, iii) nodes which are not articulation points of 
the subgraph.

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

### Output

The output is, for each pathway tested, the subnetwork, its score and its 
significance.

Need to characterize the function of this subnetwork, which is not necessarily 
limited to the pathway function. For example, how can we interpret a subnetwork 
of 10 genes inside a 1000 genes pathway? 
You have to be careful when 
interpreting the name of the pathway as a significant function beacause 
a part only of this pathway is significant.


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
