\name{signet}
% ----- \Rdversion{1.1}
\alias{signet}

\docType{package}

\title{Selection Inference in Gene Networks}

\description{
  Implements a method to perform inferences of selection in biological pathways.
  In particular it implements a simulated annealing approach that searchs for a
  high scoring subnetwork in a pathway and test for its significance.
}


\examples{
## The following example illustrates the workflow.

\dontrun{
 ## (1) The user must provide pathways and score data

 data(keggPathways);data(zScores)

 ## (2) Generate the null distribution of the subnetworks scores
 ## for subnetworks of size kmin to kmax.
 ## This may be a little long, so you can use data(nullDistExample) here.

 nullDistribution(keggPathways,zScores,iterations = 10000)

 ## (3) Apply the simulated annealing algorithm on pathways of your choice.
 ## You can provide the searchSubnet() function a graph list, or a single graph.

 searchSubnet(keggPathways[[1]],zScores,iterations = 10000)
}

}
\keyword{ package }
