if(!require("BoolNet",character.only = TRUE)) install.packages("BoolNet")
if(!require("ggplot2",character.only = TRUE)) install.packages("ggplot2")
library("BoolNet")
library("ggplot2")

starttime <- Sys.time()

every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}
update_geom_defaults("point",list(size=4))

#Running this script required 1.49h on a MacBook Pro 
#with a 2.3GHz Quad-Core processor and 16 GB of RAM

GiacoPath <- "./28nets_rules/Giacomantonio2010.txt"
FaurePath <- "./28nets_rules/Faure2006.txt"
GiacoNetwork <- loadNetwork(GiacoPath)
FaureNetwork <- loadNetwork(FaurePath)
GiacoAttr2Check <- c(0,1,0,0,1) #87.5% single state attractor
FaureAttr2Check <- c(0,1,0,0,0,1,0,1,0,0) #50% single state attractor
nrStates=10000

set.seed(123)
thetas = seq(0,180,10) #angles for biasing

#Function for converting theta bias angle to probability of activation in start state
theta2prob <- function(theta){
  theta_radians = theta * (pi/180)
  #|Rotated> := Ry(theta)|0> = cos(theta/2)|0> + sin(theta/2)|1>
  #Projection |1><1|Rotated> = sin(theta/2)|1>
  # => Probability of being expressed = amplitude^2
  probOfActivation = sin(theta_radians/2)*sin(theta_radians/2)
  return(probOfActivation)
}

#Generate nrStates start states according to given initial angles
generateStartStates <- function(initAngles, nrStates=10000){
  #need list of binary vectors of length(initAngles)
  initProbs <- sapply(initAngles, FUN=theta2prob)
  #print("initProbs:")
  #print(initProbs)
  n <- length(initAngles)
  startstates <- rep(list(NA), nrStates)
  for (s in 1:nrStates){
    state <- rep(0, n)
    randoms <- runif(n)
    for (g in 1:n){
      #overwrite entry with 1 if sample from range[0,1] is >= initProb[g]
      if (initProbs[g] >= randoms[g]){
        state[g] <- 1
      }
    }
    startstates[[s]] <- state
  }
  return(startstates)
}

#Bias initial expression of single components in Giacomantonio network
print("ANALYZING NETWORK OF GIACOMANTONIO ET AL.:")
n <- length(GiacoNetwork$genes)
thetavariationMatrix_Giaco = matrix(NA, nrow = n, ncol = length(thetas))
rownames(thetavariationMatrix_Giaco) <- GiacoNetwork$genes
rownames(thetavariationMatrix_Giaco)[5] <- "Coup-tfi"
colnames(thetavariationMatrix_Giaco) <- thetas
for (g in 1:n){
  print(paste0("Component g=",g))
  for (thetaindex in 1:length(thetas)){
    counter <- 0
    theta <- thetas[thetaindex]
    print(paste0("Angle theta=",theta))
    initangles <- rep(90, n)
    initangles[g] <- theta #Only vary bias of a single component
    #Generate list of start states
    startstates_10k = generateStartStates(initangles, nrStates=nrStates)
    
    #Follow all 10k state trajectories to their attractor individually, 
    #count how often this yields the attractor "GiacoAttr2Check", normalise by 10000
    for (s in 1:nrStates){
    pathToAttr <- getPathToAttractor(GiacoNetwork, 
                                     state=startstates_10k[[s]], 
                                     includeAttractorStates = "first")
    #Check if last row in trajectory is the desired singleton attractor
    finalState <- as.numeric(pathToAttr[nrow(pathToAttr),])
    if (identical(finalState, GiacoAttr2Check)){counter <- counter + 1}
    }
    thetavariationMatrix_Giaco[g, thetaindex] <- counter/nrStates
    print(thetavariationMatrix_Giaco)
    saveRDS(thetavariationMatrix_Giaco, file="./FigS4/thetavariationMatrix_Giaco.RDS")
  }
}

intermediatetime <- Sys.time() - starttime
print("Time after analysis of Giacomantonio network:")
print(intermediatetime)

#Bias initial expression of single components in Faure network
print("ANALYZING NETWORK OF FAURE ET AL.:")
n <- length(FaureNetwork$genes)
thetavariationMatrix_Faure = matrix(NA, nrow = n, ncol = length(thetas))
rownames(thetavariationMatrix_Faure) <- FaureNetwork$genes
colnames(thetavariationMatrix_Faure) <- thetas
for (g in 1:n){
  print(paste0("Component g=",g))
  for (thetaindex in 1:length(thetas)){
    counter <- 0
    theta <- thetas[thetaindex]
    print(paste0("Angle theta=",theta))
    initangles <- rep(90, n)
    initangles[g] <- theta #Only vary bias of a single component
    #Generate list of start states
    startstates_10k = generateStartStates(initangles, nrStates=nrStates)
    
    #Follow all 10k state trajectories to their attractor individually, 
    #count how often this yields the attractor "FaureAttr2Check", normalise by 10000
    for (s in 1:nrStates){
      pathToAttr <- getPathToAttractor(FaureNetwork, 
                                       state=startstates_10k[[s]], 
                                       includeAttractorStates = "first")
      #Check if last row in trajectory is the desired singleton attractor
      finalState <- as.numeric(pathToAttr[nrow(pathToAttr),])
      if (identical(finalState, FaureAttr2Check)){counter <- counter + 1}
    }
    thetavariationMatrix_Faure[g, thetaindex] <- counter/nrStates
    print(thetavariationMatrix_Faure)
    saveRDS(thetavariationMatrix_Faure, file="./FigS4/thetavariationMatrix_Faure.RDS")
  }
}



#Generate plot for Giacomantonio network
df_Giaco <- as.data.frame(thetavariationMatrix_Giaco)
df_Giaco["component"] <- rownames(df_Giaco)
df_Giaco.molten <- melt(df_Giaco, id.vars = "component", 
                        value.name = "AttrProb", variable.name = "theta")
df_Giaco.molten$component <- factor(df_Giaco.molten$component, 
                                    levels = c("Fgf8", "Emx2", "Pax6", "Sp8", "Coup-tfi"))
ggplot(df_Giaco.molten, aes(theta, AttrProb)) + 
  geom_point(aes(x = theta, y = AttrProb, col=component, 
                 shape=component), alpha=0.7) + 
  xlab(expression(paste("Shift in angle ", theta," [°]"))) + ylab("Probability of finding the attractor state 10010") + 
  ggtitle("c Mammalian cortical area development network \n(Classical simulation)") + 
  geom_hline(yintercept=0.875, linetype="dashed") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.85, 0.2)) +
  theme(legend.title=element_blank()) + 
  theme(legend.box.background = element_rect(colour = "black")) + 
  scale_color_manual(values = c("#e31a1c", "#cab2d6", "#33a02c", "#1f78bf", "#fdbf6f")) + 
  scale_shape_manual(values=c(16,15,25,17,8)) + 
  scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) + 
  scale_x_discrete(breaks=every_nth(n=2)) + 
  theme(axis.text = element_text(size = 25)) + 
  theme(axis.title = element_text(size = 25)) + 
  theme(title = element_text(size = 22))

ggsave(filename="./FigS4/FigS4c_Giacomantonio.pdf", 
       width = 10)

#Generate plot for Faure network
df_Faure <- as.data.frame(thetavariationMatrix_Faure)
df_Faure["component"] <- rownames(df_Faure)
df_Faure.molten <- melt(df_Faure, id.vars = "component", 
                        value.name = "AttrProb", variable.name = "theta")
df_Faure.molten$component <- factor(df_Faure.molten$component, 
                                    levels = c("CycD", "Rb", "E2F", "CycE", "CycA",
                                               "p27", "Cdc20", "Cdh1", "UbcH10", "CycB"))
ggplot(df_Faure.molten, aes(theta, AttrProb)) + 
  geom_point(aes(x = theta, y = AttrProb, col=component, 
                 shape=component), alpha=0.7) + 
  xlab(expression(paste("Shift in angle ", theta," [°]"))) + ylab("Probability of finding the attractor state 0010100010") + 
  ggtitle("d Cell cycle network \n(Classical simulation)") + 
  geom_hline(yintercept=0.5, linetype="dashed") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = c(0.85, 0.8)) +
  theme(legend.title=element_blank()) + 
  theme(legend.box.background = element_rect(colour = "black")) + 
  scale_color_manual(values = c("#e31a1c", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
                                "#a6cee3", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")) + 
  scale_shape_manual(values=c(16,15,25,17,8,10,11,12,13,14)) + 
  scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) + 
  scale_x_discrete(breaks=every_nth(n=2)) + 
  theme(axis.text = element_text(size = 25)) + 
  theme(axis.title = element_text(size = 25)) + 
  theme(title = element_text(size = 22)) + 
  theme(legend = element_text(size = 22))

ggsave(filename="./FigS4/FigS4d_Faure.pdf", 
       width = 10)

#Runtime:
runtime <- Sys.time() - starttime
print(paste0("Run time: ", runtime))
