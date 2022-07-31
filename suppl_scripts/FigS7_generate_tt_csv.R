#Create two-column parsed transition tables for all 28 analysed networks, save as csv files
if(!require("BoolNet",character.only = TRUE)) install.packages("BoolNet")
if(!require("gtools",character.only = TRUE)) install.packages("gtools")
library("BoolNet")
library("gtools")
networkpath <- "./28nets_rules/"
nets <- mixedsort(dir(networkpath, pattern = ".txt"))

#Running this script required 7.3 hours of computation time for all 28 networks
#on a MacBook Pro with a 2.3GHz Quad-Core processor and 16 GB of RAM.

starttime <- Sys.time()

#Parse transition tables into two column matrices
for (n in 1:length(nets)){
  print(paste0("Parsing transition table of net ", n, "/", length(nets)))
  network <- loadNetwork(paste0(networkpath, nets[n]))
  NrGenes <- length(network$genes)
  stateInputsInSTG <- rep(NA, 2^NrGenes)
  #tt <- readRDS(paste0("./tt_results/tt_", sub('\\.txt$', '', nets[n]), ".RDS"))
  attrs <- getAttractors(network)
  tt <- getTransitionTable(attrs)
  parsed_tt <- matrix(NA, nrow=2^NrGenes, ncol=2)
  colnames(parsed_tt) <- c("State", "Next state")
  
  for (i in 1:(2^NrGenes)){
    parsed_tt[i,1] <- paste(unlist(tt[i,1:NrGenes]), sep="", collapse="")
    parsed_tt[i,2] <- paste(unlist(tt[i,(NrGenes+1):(2*NrGenes)]), sep="", collapse="")
  }
  write.csv(parsed_tt, file = paste0("./28nets_csv_tts/parsedtt_", sub('\\.txt$', '', nets[n]), ".csv"))
}

endtime <- Sys.time()
runtime <- endtime-starttime
print("Runtime:")
print(runtime)
