# signet: Selection Inference in Gene NETworks



## Workflow

### General principle

The general idea is simple: we consider a gene list with prealably defined scores 
(e.g. a differentiation measure like the Fst) and we want to find gene networks
presenting a highest score than expected under a null hypothesis.

To do so, we will use biological pathways databases converted as gene networks 
and search in these graphs for high-scoring subnetworks.

## A walkthrough example

### Installation

There is no official release of the `signet` package at the moment. 
But you can install the development version on GitHub using the `devtools` 
package (`Rtools` must also be installed and properly configured):

```{r}
#install.packages('devtools')
devtools::install_github('CMPG/signet')
```

### Analysis

First, we generate the background distribution of the subnetworks scores 
for subnetworks of size f from 1 to 200 (the size of the biggest KEGG pathway). 
This may be a little long, so you can use `data(backgroundDist)` instead.

```{r}
backgroundDist(keggPathways,zScores,iterations = 5000)
```
Then, we apply the simulated annealing algorithm 
on pathways of your choice. Pathways must be in the `graphNEL` format. 
You can provide the `searchSubnet()` function a graph list, or a single graph.

```{r}
signetObject <- searchSubnet(pathways = keggPathways[[1]],
                             score = zScores,
                             null = backgroundDist,
                             iterations = 5000,
                             temperature = 0.995)
```

This function returns a list of N elements (corresponding to N pathways),
including a table with the whole list of genes, their scores, 
and a boolean indicating if they are found in the high-scoring subnetwork.

```{r}
testSubnet(signetObject,
           cluster = "max",
           multipleTesting = TRUE,
           threshold = 0.05)
```

The subnetwork score and the p-value are also included.

Then, you can make a correction for overlapping and multiple testing.

```{r}
results <- correctSubnet(signetObject,
                         cluster = "max",
                         multipleTesting = TRUE,
                         threshold = 0.05)
```
This will return only the subnetworks resisting to overlapping and/or 
multiple testing correction.

You can then write the results in a file in your working directory.

```{r}
writeResults(results)
```
