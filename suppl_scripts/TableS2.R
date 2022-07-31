if(!require("BoolNet",character.only = TRUE)) install.packages("BoolNet")
library("BoolNet")

file <- "./28nets_rules/Giacomantonio2010.txt"
network <- loadNetwork(file)

#Running this script required 0.02sec on a MacBook Pro 
#with a 2.3GHz Quad-Core processor and 16 GB of RAM
starttime <- Sys.time()

#List the attractors of the Giacomantonio network in its unperturbed form, 
#as well as all Pax6+Coup_tfi double perturbations
#Note that the order of bits listed here is reversed from that used in Qiskit

#Attractors of unperturbed network:
attrs <- getAttractors(network)
attrs #01001, 10110

#Attractors of OE+KO network:
perturbednetwork <- network
perturbednetwork$fixed[3] <- 1
perturbednetwork$fixed[5] <- 0
perturbedattrs <- getAttractors(perturbednetwork)
perturbedattrs #00100, 10110

#Attractors of OE+OE network:
perturbednetwork <- network
perturbednetwork$fixed[3] <- 1
perturbednetwork$fixed[5] <- 1
perturbedattrs <- getAttractors(perturbednetwork)
perturbedattrs #00101, 10111

#Attractors of KO+KO network:
perturbednetwork <- network
perturbednetwork$fixed[3] <- 0
perturbednetwork$fixed[5] <- 0
perturbedattrs <- getAttractors(perturbednetwork)
perturbedattrs #00000, 10010

#Attractors of KO+OE network:
perturbednetwork <- network
perturbednetwork$fixed[3] <- 0
perturbednetwork$fixed[5] <- 1
perturbedattrs <- getAttractors(perturbednetwork)
perturbedattrs #01001, 10011

#Runtime:
runtime <- Sys.time() - starttime
print(paste0("Run time: ", runtime))

