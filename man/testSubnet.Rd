% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testSubnet.R
\name{testSubnet}
\alias{testSubnet}
\title{Test the significance of high-scoring subnetworks found using simulated
annealing.}
\usage{
testSubnet(sigObj, null)
}
\arguments{
\item{sigObj}{A list of signet objects obtained using the
\verb{searchSubnet} function.}

\item{null}{Vector of null subnetwork scores generated using the
\code{nullDist} function.}
}
\value{
For each \verb{signet} object, a p-value is computed given the
provided emnpirical null distribution.
}
\description{
Test the significance of high-scoring subnetworks found using simulated
annealing.
}
\examples{

# Get KEGG pathways from the package graphite:
# library(graphite)
# kegg <- pathways("hsapiens", "kegg")
# kegg_human <- lapply(kegg, pathwayGraph)

data(daub13) # load the gene scores from Daub et al. (2013)

#run the search in all the pathways with 2500 iterations (default)
example <- searchSubnet(kegg_human, scores)

# generate the null distribution (here, only 5 values, but
# at least 1000 are advised)
null <- nullDist(kegg_human, scores, n = 5)
example <- testSubnet(example, null) #now, 'example' includes p-values
summary(example)
}
\keyword{annealing}
\keyword{distribution,}
\keyword{null}
\keyword{simulated}
