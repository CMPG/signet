## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ---- eval=TRUE, message=FALSE-------------------------------------------
library(signet)
library(graphite)

# pathwayDatabases() #to have a look at pathways and species available
# get the pathway list:
paths <- pathways("hsapiens", "kegg")

# convert the first 3 pathways to graphs:
kegg_human <- lapply(paths[1:3], pathwayGraph)
head(kegg_human)

## ---- eval=TRUE----------------------------------------------------------
data(daub13)
head(scores) # gene scores

## ---- eval=TRUE, message=FALSE, results="hide"---------------------------
# Run simulated annealing on the first 3 KEGG pathways:
HSS <- searchSubnet(kegg_human, scores)

## ---- echo=FALSE, eval=TRUE, message=FALSE, results="hide"---------------
null <- rnorm(1000, mean = 1.2, sd = 1)

## ---- eval=FALSE, message=FALSE------------------------------------------
#  #Generate the empirical null distribution
#  null <- nullDist(kegg_human, scores, n = 1000)

## ---- eval=TRUE----------------------------------------------------------
HSS <- testSubnet(HSS, null)

## ---- eval=TRUE----------------------------------------------------------
# Results: generate a summary table
tab <- summary(HSS)
head(tab)

# you can write the summary table as follow:
# write.table(tab,
#             file = "signet_output.tsv",
#             sep = "\t",
#             quote = FALSE,
#             row.names = FALSE)

## ---- eval=TRUE, echo=FALSE----------------------------------------------
fname <- tempfile()
writeXGMML(HSS[[1]], filename = fname)

## ---- eval=FALSE---------------------------------------------------------
#  writeXGMML(HSS[[1]], filename = "cytoscape_input.xgmml")

## ---- eval=FALSE---------------------------------------------------------
#  writeXGMML(HSS, filename = "cytoscape_input.xgmml", threshold = 0.01)

