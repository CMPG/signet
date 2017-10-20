\name{signet}
% ----- \Rdversion{1.1}
\alias{signet}

\docType{package}

\title{Selection Inference in Gene NETworks}

\description{
    Implements a simulated annealing approach to search for
    high scoring subnetworks of genes within biological pathways and to test
    for their significance.
}

\examples{
    ## The following example illustrates the workflow.

    ## (1) The user must provide pathways and score data in the same format as
    ## this example dataset:

    data(daub13)

    \dontrun {

    ## (2) Apply the simulated annealing algorithm on pathways of your choice.
    ## You can provide the searchSubnet() function a graph list,
    ## or a single graph:

    HSS <- searchSubnet(kegg_human, scores)

    ## (3) Generate the high-scores null distribution:

    null <- nullDist(kegg_human, scores)

    ## (4) Compute p-values and generate a summary table:

    HSS <- testSubnet(HSS, null)
    results <- summary(HSS)

    ## (5) Write a Cytoscape input file for visualisation:

    writeXGMML(HSS, threshold=0.01)

    }
}
\keyword{ package }
