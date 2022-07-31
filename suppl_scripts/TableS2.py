from functions import *
import numpy as np
import matplotlib.pyplot as plt
import pickle

starttime = datetime.now()

#Note that measurement is stochastic, results will differ slightly on each run
#Load dictionaries from pkl files to reproduce the plot as shown in the article
#Running this script required 0.6sec for nrshots=10000 on a MacBook Pro with a 2.3GHz Quad-Core processor and 16 GB of RAM

#Perturb two components of the Giacomantonio network
#Third compound Pax6 is perturbed with Ry(theta=3*pi/4) (partial overexpression bias), fifth compound Coup-tfi is perturbed with Ry(theta=pi/4) (partial knockout bias)
#thetatype="linear" indicates ranges of theta from [0,1] being mapped to [0, pi]

GiacoPath = "./28nets_rules/Giacomantonio2010.txt"
nrshots = 10000

#DoublePerturbedDict = allparam_exact_multiTransition_synchronous(GiacoPath, Tmax=4, nrshots=nrshots, normalizeOutput=True, sortOutput=True,
#                                                                 initActivities=[0.5,0.5,0.75,0.5,0.25], initPerturbed=[False, False, True, False, True],
#                                                                 thetatype="linear", seed_transpiler=123, seed_simulator=123)
#pickle.dump(DoublePerturbedDict, open("./TableS2/TableS2.pkl", "wb"))
DoublePerturbedDict = pickle.load(open("./TableS2/TableS2.pkl", "rb"))
print(DoublePerturbedDict)
print(outputformat(DoublePerturbedDict, digits=3))

print("Runtime = " + str(datetime.now() - starttime))


