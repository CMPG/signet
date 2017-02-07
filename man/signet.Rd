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

\dontrun{
 ## (1) The user must provide pathways and score data

 data(keggPathways);data(zScores)

 ## (2) Generate the background distribution of the subnetwork scores
 ## for different possible sizes (kmin to kmax).
 ## This may be a little long, so you can use data(nullDistExample) here.

 bkgd_dist <- backgroundDist(keggPathways,zScores)

 ## (3) Apply the simulated annealing algorithm on pathways of your choice.
 ## You can provide the searchSubnet() function a graph list, or a single graph.

 HSS <- searchSubnet(keggPathways,zScores,bkgd_dist)

 ## (4) Generate the high-scores null distribution.

 null <- searchSubnet(keggPathways,zScores,bkgd_dist)

 ## (5) Compute p-values.

 results <- summary(HSS, null)
}

}
\keyword{ package }
