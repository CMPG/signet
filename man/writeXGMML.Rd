% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/writeXGMML.R
\name{writeXGMML}
\alias{writeXGMML}
\title{Write Cytoscape input file}
\usage{
writeXGMML(sigObj, filename = "signet_output.xgmml", threshold = 0.01)
}
\arguments{
\item{sigObj}{A signet or signetList object.}

\item{filename}{The desired file name. Default is "signet_output.xgmml".}

\item{threshold}{Significance threshold (default: 0.01). If a signetList is
provided, all subnetworks with a p-value below this threshold will be merged
and represented.}
}
\value{
Writes an XGMML file in the working directory. If a single pathway
(signet object) is provided, the whole pathway is represented and nodes
belonging to the highest-scoring subnetwork (HSS) are highlighted in red. If
a list of pathways (signetList) is provided, all subnetworks with a p-value
below a given threshold (default: 0.01) are merged and represented. Note that
in this case, only the nodes belonging to HSS are kept for representation.
}
\description{
This function allows to write an XGMML file to represent the results in
Cytoscape.
}
\examples{
# Get KEGG pathways from the package graphite:
# library(graphite)
# kegg <- pathways("hsapiens", "kegg")
# kegg_human <- lapply(kegg, pathwayGraph)

data(daub13) # load the gene scores from Daub et al. (2013)

#run the search in all the pathways with 2500 iterations (default)
example <- searchSubnet(kegg_human, scores)

#write Cytoscape input file for the first pathway:
#writeXGMML(example[[1]], filename=tempfile())
}
\keyword{Cytoscape,}
\keyword{network,}
\keyword{plot,}
\keyword{subnetwork,}
\keyword{visualization}
