## command line usage
### Rscript iGraphNetworkRunner.R [outcome] [predictors] [zed] [method] [cutoff] [perturb] [force] [cores] [output] [cluster] [comms]
### outcome = 2 columns ... ID and PHENO
### predictors = data used to predict the pheno, must be numeric / continuous for best results
### zed = "yes" or "no", do you want your data Z transformed after filtering to make everything on the same scale
### method = "effect" or "p"
### cutoff = decimal from 0-1, representing the minimum p-value for inclusion of a feature or quantile of effect estimate you want included int he filtering (ie in this case 0.25 will include the top 25% |beta|)
### perturb = "no" or file containing list of feature names in a single column with no header to be used to make exclusions from data in analysis to simulate knock outs
### force = "no" or file containing list of feature names in a single column with no header to be used to force inclusion into network estimates
### cores = 1 - 24, how many cores you want this to run on
### output = text string to uniquely identify this run
### cluster = "corr" or "euc" for correlation or euclidian
### comms = "edge" or "spin" or "walk", denoting edge betweens, Potts spinglass method or random walk to identify communities

## set the environment and options

### packages needed
library("data.table")
library("parallel")
library("igraph")

### arguments from command lines
args <- commandArgs()
print(args)
outcomeFile <- args[6] # outcome file
predictorFile <- args[7] # predictors file
zed <- args[8] # "yes" or "no" to the Z transform
method <- args[9] # "effect" or "p" for feature filtering
cutoff <- as.numeric(args[10]) # 0-1
perturb <- args[11] # perturb predictors file to exclude
force <- args[12] # force predictors file to include
cores <- as.numeric(args[13]) # 1-24
output <- args[14]
clusterType <- args[15]
comms <- args[16]

# ### for testing
# setwd("~/Desktop/scratch/igraphNetworks/test/")
# outcomeFile <- "ppmiForNetwork_pheno.tab" # outcome file
# predictorFile <- "ppmiForNetwork_geneExp_symbolNames.tab" # predictors file
# zed <- "no" # "yes" or "no" to the Z transform
# method <- "p" # "effect" or "p" for feature filtering
# cutoff <- 0.01 # 0-1
# perturb <- "no" # perturb predictors file to exclude
# force <- "no" # force predictors file to include
# cores <- 2 # 1-24
# output <- "tuesday"
# clusterType <- "euc"
# comms <- "walk"

### test CMD
# cd ~/Desktop/scratch/igraphNetworks/test/
# Rscript iGraphNetworkRunner.R ppmiForNetwork_pheno.tab ppmiForNetwork_geneExp_symbolNames.tab no p 0.001 no no 2 testMonay1030 euc spin

## load up the data and start making thresholds for filtering

### data incoming
outcomeData <- fread(file = outcomeFile, header = T)
predictorData <- fread(file = predictorFile, header = T)
preFilteredData <- merge(outcomeData, predictorData, by = "ID")

## here is where the filtering starts

### build mcapply loop to generate quick summary stats for each predictor
predictorList <- names(predictorData)[2:length(names(predictorData))]
runPreFiltering <- function(i)
{
  thisPredictor <- as.character(predictorList[i])
  beta <- NA
  se <- NA
  p <- NA
  try({testingThePredictor <- lm(PHENO ~ preFilteredData[[thisPredictor]] , data = preFilteredData)
  beta <- summary(testingThePredictor)$coefficients["preFilteredData[[thisPredictor]]","Estimate"]
  se <- summary(testingThePredictor)$coefficients["preFilteredData[[thisPredictor]]","Std. Error"]
  p <- summary(testingThePredictor)$coefficients["preFilteredData[[thisPredictor]]","Pr(>|t|)"]})
  c(thisPredictor,beta,se,p)
}

### runs the code above then this step here exports your preFiltered association summary stats between your predictors and your outcome
nRuns <- length(predictorList)
outPut <- mclapply(1:nRuns, runPreFiltering, mc.cores = cores)
newResult <- do.call(rbind,outPut)
outfile <- paste("preFilteringSummary_",output,".tab", sep = "")
write.table(newResult, file = outfile, sep = "\t", row.names = F, quote = F, col.names = c("PREDICTOR","BETA","SE","P")) 

## now let's filter

### read back in the preFilteringSummary data
tempFilterFile <- fread(file = outfile, header = T)

### now drop anything with missing summary stats
### then generate the absolute value of effect estimates and rank by that
filterFile <- subset(tempFilterFile, !is.na(tempFilterFile$P) & !is.na(tempFilterFile$BETA))
filterFile$absBeta <- abs(filterFile$BETA)
filterFile$quantileRankedAbsBeta <- rank(filterFile$absBeta)/length(filterFile$absBeta)

### here we filter by either effect estiamte ranked quantile or P
### if excludedFromAnalysis == 1, this predictor will not be included
if(method == "effect")
{
  filterFile$excludedFromAnalysis <- ifelse(filterFile$quantileRankedAbsBeta > cutoff, 0, 1)
}
if(method == "p")
{
  filterFile$excludedFromAnalysis <- ifelse(filterFile$P < cutoff, 0, 1)
}

### sometime you want to force certain genes in or out, we accomodate those here
if(force != "no")
{
  toForce <- fread(file = force, header = F)
  for(i in 1:length(toForce$V1))
  {
    try({predictorToForce <- as.character(toForce$V1[i])
    filterFile$excludedFromAnalysis[filterFile$PREDICTOR == predictorToForce] <- 0})
  }
}
if(perturb != "no")
{
  toPerturb <- fread(file = perturb, header = F)
  for(i in 1:length(toPerturb$V1))
  {
    try({predictorToPerturb <- as.character(toPerturb$V1[i])
    filterFile$excludedFromAnalysis[filterFile$PREDICTOR == predictorToPerturb] <- 1})
  }
}

### now export a summary of filtered variants
filteredFile <- subset(filterFile, excludedFromAnalysis == 0)
outfile <- paste("filteredPredictorList_",output,".tab", sep = "")
write.table(filteredFile, file = outfile, sep = "\t", row.names = F, quote = F) 

## this is where the formatting for iGraph begins

### make reduced predictor dataset
keepPredictorsList <- filteredFile$PREDICTOR
predictorsToPull <- c("ID",keepPredictorsList)
reducedData <- preFilteredData[, ..predictorsToPull]
outfileReduced <- paste("reducedData_",output,".tab", sep = "")
write.table(reducedData, file = outfileReduced, sep = "\t", row.names = F, quote = F) 

### if zed was specified then also make a z transformed version
### zed option which is useful in heterogenous data sources and  we build all new data
### also, if you didn't know "Zed's dead baby"
### sets anything by default to zero unless the data gets fit
if(zed == "yes")
{
  toZed <- fread(outfileReduced, header = T)
  zeding <- matrix(nrow = length(toZed$ID), ncol = length(names(toZed)), NA)
  for(i in 1:length(names(toZed)))
  {
    if(i == 1)
    {
      zeding[,i] <- toZed$ID
    }
    if(i > 1)
    {
      zeding[,i] <- (toZed[[i]] - mean(toZed[[i]], na.rm = T))/sd(toZed[[i]], na.rm = T)
    }
  }
  outfileZed <- paste("reducedData_Ztransformed_",output,".tab", sep = "")
  write.table(zeding, file = outfileZed, sep = "\t", row.names = F, quote = F, col.names = names(toZed)) 
}

## now after over 160 lines of code, we finally build the damn networks

### read in the reduced data
if(zed == "yes")
{
  data <-  fread(outfileZed, header = T)
}
if(zed == "no")
{
  data <- fread(outfileReduced, header = T)
}

## format the data
### flip it and make the row names the predictors and the column names the samples
transposed <- t(data)
outfileToCluster <- paste("reducedData_toCluster_",output,".tab", sep = "")
write.table(transposed, file = outfileToCluster, sep = "\t", row.names = T, quote = F, col.names = F) 

### now read it back in
toCluster <- fread(outfileToCluster, header = T)
row.names(toCluster) <- toCluster$ID
toCluster$ID <- NULL


## make clusters based on clusterType
if(clusterType == "euc")
{
  g <- graph.adjacency(
    as.matrix(dist(toCluster, method="euclidean")),
    mode="undirected",
    weighted=TRUE,
    diag=FALSE)
}

if(clusterType == "corr")
{
  g <- graph.adjacency(
    as.matrix(as.dist(cor(t(toCluster), method="pearson"))),
    mode="undirected",
    weighted=TRUE,
    diag=FALSE)
}

## simplify clustering and clean up

### simplify
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)

### Convert edge weights to absolute values
E(g)$weight <- abs(E(g)$weight)

### Remove edges below absolute weight 0.8
g <- delete_edges(g, E(g)[which(E(g)$weight<0.8)])

### Assign names to the graph vertices, ie predictor rownames
V(g)$name <- rownames(toCluster)

### Convert the graph adjacency object into a minimum spanning tree
mst.prim <- mst(g, algorithm="prim")

### make communiites
if(comms == "edge")
{
  mst.communities <- edge.betweenness.community(mst.prim, weights=NULL, directed=FALSE)
}
if(comms == "spin")
{
  mst.communities <- spinglass.community(mst.prim, weights = NULL)
}
if(comms == "walk")
{
  mst.communities <- walktrap.community(mst.prim, weights = NULL)
}

mst.clustering <- make_clusters(mst.prim, membership=mst.communities$membership)

### Check the vertex degree, i.e., number of connections to each vertex
degree(mst.prim)

### Output information for each community, including vertex-to-community assignments and modularity
commSummary <- data.frame(
  mst.communities$names,
  mst.communities$membership,
  mst.communities$modularity
)
colnames(commSummary) <- c("Feature", "Community", "Modularity")
commSummaryOut <- commSummary[order(commSummary$Community,commSummary$Modularity),]
write.table(commSummaryOut, file = paste("networkClusters_",output,".tab", sep =""), sep = "\t", row.names = F, quote = F)


## plotting related stuff below

### overplotting here is common when 100s of predictors are included (this is a work in progress)

### community plot
pdf(file = paste("clusters_",output,".pdf", sep = ""))
plot(
  mst.clustering, mst.prim,
  vertex.label.cex=0.25
  )

dev.off()

## turn it off

q("no")
