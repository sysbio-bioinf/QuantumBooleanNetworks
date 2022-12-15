#Generate boxplots for s steps backwards up to maxstepsback across all attrs in all 28 networks
# Plot number of predecessors M/2^netsize
#   Get mean/median/SD/IQR for each step across nets
if(!require("BoolNet",character.only = TRUE)) install.packages("BoolNet")
if(!require("gtools",character.only = TRUE)) install.packages("gtools")
if(!require("ggplot2",character.only = TRUE)) install.packages("ggplot2")
if(!require("scales",character.only = TRUE)) install.packages("scales")
library("BoolNet")
library("gtools")
library("ggplot2")
library("scales")

#Running this script with the provided csv files required 14.4 minutes 
#on a MacBook Pro with a 2.3GHz Quad-Core processor and 16 GB of RAM.

#The generation of the csv files with the transition tables is specified in the generate_tt_csv.R script.
#These files need to be generated first in order to create the plots specified below.
#Note that this required 7.3 hours on the same MacBook Pro and the resulting set 
#of 28 transition tables has a total size of 2.33GB.

starttime <- Sys.time()

conv2adjmat <- function(sbmlnet, inputcorrected = FALSE){
  "Converts an SBML object to adjacency matrix. 
  If inputcorrected = T, all input edges are set to zero. 
  A vertex v is said to be an input if it is only regulated by itself, 
  meaning the sum of column v in the adjacency matrix is one, with the only non-zero entry 
  being at position [v,v]."
  adjmat <- sapply(sbmlnet$interactions, function(gene) {v <- rep(0,length(sbmlnet$genes)); 
  v[gene$input] <- 1; return(v)})
  if (inputcorrected == TRUE){
    for (d in 1:dim(adjmat)[1]){
      if (adjmat[d,d] >= 1){adjmat[d,d] <- 1}
    }
    for (col in 1:dim(adjmat)[2]){
      if (sum(adjmat[,col]) == 1 & adjmat[col,col] == 1){
        adjmat[col,col] <-0
      }
    }
  }
  return(adjmat)
}

maxstepsback <- 10
networkpath <- "./28nets_rules/"
csvpath <- "./28nets_csv_tts/"
nets <- mixedsort(dir(networkpath, pattern = ".txt"))
csvs <- mixedsort(dir(csvpath, pattern = ".csv"))
averagePredecessorPercentagesAcrossNets <- matrix(NA, nrow = length(nets), ncol = maxstepsback)

networksizes <- networkedges <- rep(NA, length(nets))
for (n in 1:length(nets)){
  network <- loadNetwork(paste0(networkpath, nets[n]))
  networksizes[n] <- length(network$genes)
  adjmat <- conv2adjmat(network)
  networkedges[n] <- sum(adjmat)
}
names(networksizes) <- names(networkedges) <- sub('\\.txt$', '', nets) 

for (n in 1:length(nets)){
  print(paste0("NETWORK n=", n, "/", length(nets)))
  statespacesize <- rep(2^networksizes[n], maxstepsback)
  parsed_tt <- read.csv(file=paste0(csvpath, csvs[n]), sep = ",", header = TRUE, colClasses=rep("character",3))[,2:3]
  
  #Get mean across all rows, i.e. over all possible markedState attractor starting points
  #Generate predecessorMatrix from csv for this network
  
  network <- loadNetwork(paste0(networkpath, nets[n]))
  attrs <- getAttractors(network, method="sat.exhaustive")
  attrstatesvector <- c()
  for (attr in attrs$attractors){
    invStates <- attr$involvedStates #Will always have only one row since all networks have netsize < 32
    for (state in invStates){
      attrbitstring <- paste(as.character(BoolNet:::dec2bin(state, 32)[1:networksizes[n]]), collapse="")
      attrstatesvector <- append(attrstatesvector, attrbitstring)
    }
  }
  #print("All attractor states in this network:")
  #print(attrstatesvector) #Vector with bitstrings of all attractor states
  nrattractorstates <- length(attrstatesvector)
  predecessorMatrix <- matrix(NA, nrow = nrattractorstates, ncol = maxstepsback) #matrix to be filled and saved
  
  for (i in 1:nrattractorstates){
    print(paste0("Attractor state ", i, " out of ", nrattractorstates))
    attrstate <- attrstatesvector[i]
    solutionstates <- attrstate
    for (s in 1:maxstepsback){
      print(paste0("Step back = ", s, " out of ", maxstepsback))
      successorIndices <- which(parsed_tt[,2] %in% solutionstates)
      solutionstates <- parsed_tt[successorIndices,1]
      predecessorMatrix[i,s] <- length(solutionstates) #row=attr, col=nr steps back
    }
  }
  print(predecessorMatrix)
  #Columns of predecessorMatrix will come closer to summing up to N=2^n 
  #as the number of backward steps increases and GoE states are reached
  saveRDS(predecessorMatrix, 
          file=paste0("predecessorMatrices/net_", 
                      n, "_maxstepsback_", maxstepsback, "_predecessorMatrix.RDS"))
  
  meanOverAttrs <- colMeans(predecessorMatrix)
  meanOverAttrsPercentage <- meanOverAttrs/statespacesize
  averagePredecessorPercentagesAcrossNets[n,] <- meanOverAttrsPercentage
}
print(averagePredecessorPercentagesAcrossNets)

#Save data:
#saveRDS(averagePredecessorPercentagesAcrossNets, 
#        file="predecessorMatrices/averagePredecessorPercentagesAcrossNets.RDS")

#Load data:
averagePredecessorPercentagesAcrossNets <- readRDS(file="predecessorMatrices/averagePredecessorPercentagesAcrossNets.RDS")

df <- ylogdf <- as.data.frame(matrix(NA, nrow = 0, ncol = 3))
colnames(df) <- colnames(ylogdf) <- c("avgfractionacrossattrs", "net", "step")
for (n in 1:dim(averagePredecessorPercentagesAcrossNets)[1]){
  for (s in 1:dim(averagePredecessorPercentagesAcrossNets)[2]){
    df[nrow(df) + 1,] = c(averagePredecessorPercentagesAcrossNets[n,s], n, s)
    ylogdf[nrow(ylogdf) + 1,] = c(log(averagePredecessorPercentagesAcrossNets[n,s]), n, s)
  }
}

pdf(file="FigS7/FigS7_linear.pdf")
ggplot(df, aes(x=step, y=avgfractionacrossattrs, group=step)) + geom_boxplot(color="black", fill="white") + 
  theme_bw(base_size = 22) + theme(panel.background = element_blank()) + 
  xlab("Number of inverted transitions") + 
  ylab("Average M/N across attractor states\n for all networks") + + 
  scale_y_continuous(breaks = seq(0, max(df$avgfractionacrossattrs), by = 0.05)) + 
  scale_x_continuous(breaks = seq(min(df$step), max(df$step), by = 1))
dev.off()


scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
pdf(file="FigS7/FigS7_logaxis.pdf")
ggplot(df, aes(x=step, y=avgfractionacrossattrs, group=step)) + geom_boxplot(color="black", fill="white") + 
  theme_bw(base_size = 22) + theme(panel.background = element_blank()) + 
  xlab("Number of inverted transitions") + 
  ylab("Average M/N across attractor states\n for all networks") + 
  scale_x_continuous(breaks = seq(min(df$step), max(df$step), by = 1)) + 
  scale_y_log10(label=scientific_10)
dev.off()

print("Median and IQR values for every number of inverted transitions:")
median(averagePredecessorPercentagesAcrossNets[,1])
apply(averagePredecessorPercentagesAcrossNets, 2, median)
apply(averagePredecessorPercentagesAcrossNets, 2, IQR)

print("Avg+SD of network nodes/edges:")
round(mean(networksizes),1)   #15.5
round(sd(networksizes),1)   #5.4
round(mean(networkedges),1) #38.6
round(sd(networkedges),1)   #17.1

#Get correlations for average M/N by network for s=1 with networksizes
print("Correlations between network size and M/N after one inverted step:")
print(round(cor(networksizes, averagePredecessorPercentagesAcrossNets[,1]),3)) #-0.714
#Larger networks tend to have a smaller fraction M/N after s=1

endtime <- Sys.time()
runtime <- endtime - starttime
print(paste0("Runtime = ", runtime))

#"Median and IQR values for every number of inverted transitions:"
#MEDIAN:
#[1] 0.003021967 0.015349153 0.043443468 0.059497070 0.075668570 0.082865212 0.096875000
#[8] 0.099414063 0.100000000 0.100000000
#IQR:
#[1] 0.02267585 0.05563202 0.08830041 0.12030986 0.11523292 0.12286795 0.12797964 0.12700181
#[9] 0.12575059 0.12395824

