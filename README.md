# signet: Selection Inference in Gene NETworks

We have developed a method to detect selection in biological pathways. 
The general idea is to search for gene subnetworks within biological pathways
that present unusual features. This search is a typical combinatorial 
optimization problem that can be solved using a heuristic approach like 
simulated annealing. We implemented such an algorithm to search for high-scoring
subnetworks of genes within biological pathways. We also developed a 
significance testing procedure that explicitly takes into account the 
optimization process involved in this search.

This workflow is implemented in the `signet` package presented here.

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
### Input

`signet` takes as input a data frame of gene scores. the first column
The other input is a list a graph representing biological pathways. We advise 
to use the package `graphite` to generate this list:

```{r}
library(graphite)
pathwayDatabases() #check pathways and species available
paths <- pathways("hsapiens", "kegg") #get the pathway list
graphs <- lapply(paths, pathwayGraph) #convert pathways to graphs
```
Note that gene identifiers must be the same between the gene scores data frame
and the pathway list (e.g. entrez). `graphite` provides a function to convert
gene identifiers.

A example dataset from Daub et al. (2013) as well as human KEGG pathways are 
provided:

```{r}
data(daub13)
```
### Workflow

We first have to generate a "background distribution" of subnetworks scores
to be able to compare subnetworks of different sizes with simulated annealing.
To do so, we sample random subnetworks and compute the distribution of 
their scores for different possible sizes:
```{r}
#Generate background distribution:
bkgd_dist <- backgroundDist(kegg_human,scores)
```
This distribution is then used to search for high-scoring subnetworks within the
provided biological pathways:

```{r}
#Run simulated annealing on the first 10 pathways:
HSS <- searchSubnet(kegg_human[1:10],
                    scores,
                    bkgd_dist,
                    iterations = 5000)
```
This function returns, for each pathway, the highest-scoring subnetwork found,
its score and other information about the simulated annealing run. You can have
a look at the network and the optimization process as follow:

```{r}
#Inspect a pathway:
plot(HSS[[1]])
```
Then, to test the significance of the high-scoring subnetworks, you have to 
generate a null distribution of high-scores:

```{r}
#Generate the null distribution
null <- nullDist(kegg_human,scores,1000,bkgd_dist)
```
Note that the `null` object is a simple vector of null high-scores. Therefore,
you can run other iterations afterwards and concatenate the output with the
previous vector to compute more precise p-values.

This distribution is finally used to compute p-values and update the 
`signet` object:

```{r}
#Compute p-values
HSS <- testSubnet(HSS,null)
```
### Interpretation of the results

When p-values have been computed, you can generate a summary table
(one row per pathway):
```{r}
#Results: generate a summary table
tab <- summary(HSS)

```
Note that searching for high-scoring subnetworks and generating the null 
distribution usually takes a few hours. However, these steps are easy
to parallelize on a cluster as different iterations are independent from each 
other. Depending on the number of processors available, the computation time
can be reduced to a few minutes.


