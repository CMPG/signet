# signet: Selection Inference in Gene NETworks

We have developed a method to detect selection in biological pathways. 
The general idea is to search for subsets of interacting genes within gene 
networks (e.g. biological pathways) that present unusual features. This search 
is a typical combinatorial optimization problem that can be solved using a 
heuristic approach like simulated annealing. We implemented such an algorithm 
to search for high-scoring subnetworks of genes in biological pathways data. 
We also developed a significance test procedure that explicitly takes into 
account the optimization process involved in this search.

This workflow is implemented in the `signet` package presented here. It is an 
extension of genome scans for selection to gene networks. Using a statistic to 
measure selection, the idea is to find in gene networks high scoring subnetworks 
of genes within biological pathways.

## Workflow

### General principle

The general idea is simple: we consider a gene list with prealably defined scores 
(e.g. a differentiation measure like the Fst) and we want to find gene networks
presenting a score higher than expected under the null hypothesis.

To do so, we will use biological pathways databases converted as gene networks 
and search in these graphs for high-scoring subnetworks.

<p align="center"><img src="misc/workflow.png" width="400"></p>
<p align="center">Figure 1: Workflow.</p>

### Input preparation

#### Biological networks conversion to graphs

In order to get information about relationships between genes,
we advise to use the Bioconductor package `graphite`, which implements a 
procedure to convert biological pathways data from different databases to 
graphs (for more details, see Sales et al., 2012). These databases have emerged 
as references in systems biology: KEGG, Reactome, BioCarta, NCI, Panther 
and HumanCyc.

#### Gene scores


### Search algorithm

The method implemented is based on Ideker et al. (2002) heuristics, 
but several improvements are considered. You can see below an animation 
representing a run of the simulated annealing algorithm used in the package.

<p align="center"><img src="misc/anim_50fps.gif"></p>
<p align="center">Figure 2: Simulated annealing run in a gene network.</p>

As you can see, as we add or remove new genes in the active subnetwork (in red),
and the subnetwork score is maximized as we iterate.

#### Background distribution

To apply the algorithm to a list of pathways, you have to first generate the 
background distribution of your summary statistic.

The input is the list of all the pathways you want to analyze, and a list of 
scores for each gene. The background distribution or the statistic will be 
generated for each possible subnetwork size (kmin to kmax).

A pathway is randomly sampled (the sampling probability being conditioned by 
the number of genes in the pathway) and a connected subnetwork of size k is 
randomly picked (one gene is sampled in the pathway, then k genes in the 
boundary). The score of the subnetwork is then computed. Here, the score is 
simply the average Fst over all genes in the subnetwork.

This is done N times for each k, to get the background distribution of 
subnetworks scores, we just keep the mean and standard error of this 
distribution.

#### Algorithm for high-scoring subnetwork search

We consider that genes can yield two states: active or inactive.

1. Generate a random solution, i.e. pick a random subnetwork of arbitrary size k
in the gene network.

2. Calculate its score using a scoring function. Here, the 
score s is the average FST, standardized by the background distribution for
a network of size k.

3. Generate a random neighboring solution, i.e. a new subnetwork.
To do so, we randomly pick a new gene and toggle its state (i.e. we activate
or inactivate it). The selected gene is not an articulation point of 
the subgraph, i.e. its removal doesn’t disconnect the active subgraph. The 
selected gene is randomly picked from the following nodes: i) nodes in 
the boundary; ii) leaves, iii) nodes which are not articulation points of 
the subgraph.

4. Calculate the new subnetwork's score

5. Compare them:
If snew < sold: move to the new solution
If snew > sold: maybe move to the new solution (acceptance probability)

6. Repeat steps 3-5 above until an acceptable solution is found or you reach 
some maximum number of iterations.

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

### Output

The output is, for each pathway tested, the high-scoring subnetwork, 
its score and its level of significance.

## A walkthrough example

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
