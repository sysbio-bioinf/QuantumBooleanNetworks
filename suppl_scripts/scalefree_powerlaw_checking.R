#Testing networks for scale-free topology using the poweRlaw package
#Networks are retained given p-values above a threshold of p=0.1,
#meaning that the power law is a plausible hypothesis for the data, 
#as described by Clauset et al. (2009), Power-Law Distributions in Empirical Data

#Running this script required 1.06mins on a MacBook Pro 
#with a 2.3GHz Quad-Core processor and 16 GB of RAM

if(!require("poweRlaw",character.only = TRUE)) install.packages("poweRlaw")
if(!require("gtools",character.only = TRUE)) install.packages("gtools")
library("poweRlaw",character.only = TRUE)
library("gtools")

starttime <- Sys.time()

readNetwork <- function(nwFileName, nwName) {
  inFile <- paste(nwFileName,sep="")
  return(read.csv(inFile, stringsAsFactors=FALSE))
}

simplify <- function(genes, rules){
  srules <- list()
  for (r  in 1:length(rules)) {
    rule <- gsub(" ", "", rules[r])
    rule <- gsub("&", ",",rule, fixed=TRUE)
    rule <- gsub("|", ",",rule, fixed=TRUE)
    rule <- gsub("(", "", rule, fixed=TRUE)
    rule <- gsub(")", "", rule, fixed=TRUE)
    rule <- gsub("!","", rule, fixed=TRUE)
    srules[r] <- list(which(genes %in% strsplit(rule,",")[[1]]))
  }
  return (srules)
}

countDegree <- function(network, name, noOfSims = 500, noOfThreads=8 ) {
  ccNw <- readNetwork(network, name)
  ccGenes <- ccNw[,1]
  ccRules <- ccNw[,2]
  sr <- simplify(ccGenes,ccRules)
  degs <- matrix(0,ncol=6, nrow=length(ccGenes), dimnames=list(ccGenes,list("out","inp","loop","total","Z (out)", "Z (total)")))
  degs[,1] <- sapply(1:length(ccGenes),function(g) { length(sr[[g]]) })
  degs[,2] <- c(tabulate(unlist(sr)), rep(0,length(ccGenes)-max(unlist(sr))))
  degs[,3] <- sapply(1:length(ccGenes), function(g) { if (g %in% sr [[g]]) {return(1)} else {return(0)} } )
  degs[,1] <- degs[,1] - degs[,3]
  degs[,2] <- degs[,2] - degs[,3]
  degs[,4] <-  degs[,1] + degs[,2] + degs[,3]
  degMean <- mean(degs[,1])
  degSd <- sd(degs[,1])
  degs[,5] <- round((degs[,1] - degMean)/degSd,2)
  degMean <- mean(degs[,4])
  degSd <- sd(degs[,4])
  degs[,6] <- round((degs[,4] - degMean)/degSd,2)
  degtable <- degs
  degtable[degtable==0] <- NA
  m <- degs[,4]
  m <- m[m>0]
  m_pl <- displ$new(m)
  est <- estimate_xmin(m_pl)
  if(is.na(est$xmin)) {
    print("failed to set model parameters")
    pvalue <- "failed to set model parameters"
  } else {
    m_pl$setXmin(est)
    bt_pl <- bootstrap_p(m_pl, no_of_sims=noOfSims, threads=noOfThreads)
    pvalue <- bt_pl$p
    print(paste0("poweRlaw bootstrap_p p-value = ", pvalue))
  }
  return(list(bootstrap=bt_pl, degrees=degs))
}


networkpath <- "./28nets_rules/"
nets <- mixedsort(dir(networkpath, pattern = ".txt"))
for (n in 1:length(nets)){
  print(paste0("Network = ", nets[n]))
  network <- paste0(networkpath, nets[n])
  Zscore <- countDegree(network,network, noOfSims=100, noOfThreads=4)
}
#All networks yielded a p-value > 0.1 -> power law is plausible
runtime <- Sys.time() - starttime
print(paste0("Run time: ", runtime))
