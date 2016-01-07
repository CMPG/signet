### Example

# library(roxygen2)
# library(devtools)
# document()

library(signet)
data(keggPathways)
data(zScores)

bkgd <- backgroundDist(pathwaysList = keggPathways,
                       scores = zScores,
                       kmax = 10)
# data(bkgd)


signetObject <- searchSubnet(pathway = keggPathways,
                             scores = zScores,
                             nullDist = bkgd)


results <- correctOverlap(signetObject,
                          cluster = "max",
                          multipleTesting = TRUE,
                          threshold = 0.05)

writeResults(results)
