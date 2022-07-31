from functions import *
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import copy
from datetime import datetime
import pickle

starttime = datetime.now()
pickledump = True #Should results of this run be saved (True), else load previously stored results (False)

#Load dictionaries from pkl files to reproduce the plot as shown in the article
#Running this script required 10min18sec for shotnr=10000 on a MacBook Pro with a 2.3GHz Quad-Core processor and 16 GB of RAM

FaurePath = "./28nets_rules/Faure2006.txt"
shotnr=10000


#RUN STATE TRANSITION CIRCUITS
#T=1
FaureExactTransition1 = allparam_exact_multiTransition_synchronous(rulestxt=FaurePath, Tmax=1, nrshots=shotnr, normalizeOutput=True, seed_transpiler=123, seed_simulator=123)
print(FaureExactTransition1)
print(len(FaureExactTransition1.keys()))
#T=9
FaureReinitTransition9 = allparam_multiTransition_synchronous(rulestxt=FaurePath, Tmax=9, nrshots=shotnr, normalizeOutput=True, seed_transpiler=123, seed_simulator=123)
print(FaureReinitTransition9)
print(len(FaureReinitTransition9.keys()))
#45.92% in single state attractor, 54.08% across the remaining 7 states of the cyclic attractor

if pickledump:
    pickle.dump(FaureExactTransition1, open("./FigS4/Faure_T1dictionary.pkl", "wb"))
    pickle.dump(FaureReinitTransition9, open("./FigS4/Faure_T9dictionary.pkl", "wb"))
else:
    FaureExactTransition1 = pickle.load(open("./FigS4/Faure_T1dictionary.pkl", "rb"))
    FaureReinitTransition9 = pickle.load(open("./FigS4/Faure_T9dictionary.pkl", "rb"))


print("Begin calculation for initial theta variation in Faure network... ")
FullSyncCirc = synthesizeFullNetworkUpdateCircuit(rulestxt=FaurePath, update="synchronous")
n = FullSyncCirc.num_qubits//2
Tmax_Faure = 9
allthetas = list(range(0,190, 10)) #0,...,180
thetainitvec = list(np.repeat(90, n))

resultmatrix = np.zeros((n, len(allthetas)), dtype=float)

for g in range(n):
    print("Gene " + str(g))
    for thetanr in range(len(allthetas)):
        thetashift = allthetas[thetanr]
        print("thetashift " + str(thetashift))
        # Change to only vary thetainitvec in gene g
        thetainitvec_g = copy.deepcopy(thetainitvec)
        thetainitvec_g[g] = thetashift

        Transition = allparam_multiTransition_synchronous(rulestxt=FaurePath, Tmax=Tmax_Faure,
                                                         nrshots=shotnr, initActivities=thetainitvec_g,
                                                         thetatype="angle", normalizeOutput=True,
                                                         seed_transpiler=123, seed_simulator=123)

        if "0010100010" in Transition.keys():
            resultmatrix[g,thetanr] = Transition["0010100010"]
        else:
            print("Attractor not measured!")
            resultmatrix[g,thetanr] = 0

        #Save results
        if pickledump:
            np.save(file="./FigS4/FaureShiftByTheta.npy", arr=resultmatrix)
resultmatrixFaure =  resultmatrix

#Load previous results
if pickledump == False:
    resultmatrixFaure = np.load(file="./FigS4/FaureShiftByTheta.npy")


#Plot variation in Giacomantonio network
#Load theta variation data calculation for main figure 2b
resultmatrix_Giaco = np.load(file="./FigS4/GiacomantonioShiftByTheta.npy")

print("Begin plotting Giacomantonio network... ")
plt.scatter(allthetas, resultmatrix_Giaco[0,:], c="#e31a1c", label="Fgf8", edgecolors='none', alpha=0.5, marker="o")
plt.scatter(allthetas, resultmatrix_Giaco[1,:], c="#cab2d6", label="Emx2", edgecolors='none', alpha=0.75, marker="s")
plt.scatter(allthetas, resultmatrix_Giaco[2,:], c="#33a02c", label="Pax6", edgecolors='none', alpha=0.5, marker="v")
plt.scatter(allthetas, resultmatrix_Giaco[3,:], c="#1f78b4", label="Sp8", edgecolors='none', alpha=0.75, marker="^")
plt.scatter(allthetas, resultmatrix_Giaco[4,:], c="#fdbf6f", label="Coup-tfi", edgecolors='none', alpha=1, marker="*")
plt.legend(loc='lower right')
plt.grid(False)
plt.ylim((-0.05,1.05))
plt.xlabel(r'Shift in angle $\theta$ [°]', fontsize=12)
plt.ylabel("Probability of measuring the attractor state 10010", fontsize=12)
#plt.title(r'Deviation from unbiased single state basin size by varying angle $\theta$')
#plt.title(r'$\bf{a}$ Mammalian cortical area development network')
plt.title('$\\bf{a}$ Mammalian cortical area development network \n(Quantum simulation)')
ax = plt.gca()
plt.xticks(allthetas)
labelnr=0
for label in ax.xaxis.get_ticklabels():
    if labelnr % 2 != 0:
        label.set_visible(False)
    labelnr += 1

plt.axhline(y=0.875, color="black", xmin=0, xmax=180, linestyle="dashed")
plt.savefig('./FigS4/FigS4a_Giacomantonio.pdf', dpi=350)
plt.clf()

#Plot variation in Faure network
print("Begin plotting Faure network... ")
plt.scatter(allthetas, resultmatrixFaure[0,:], c="#e31a1c", label="CycD", edgecolors='none', alpha=0.5, marker="o")
plt.scatter(allthetas, resultmatrixFaure[1,:], c="#1f78b4", label="Rb", edgecolors='none', alpha=0.5, marker="s")
plt.scatter(allthetas, resultmatrixFaure[2,:], c="#b2df8a", label="E2F", edgecolors='none', alpha=0.5, marker="v")
plt.scatter(allthetas, resultmatrixFaure[3,:], c="#33a02c", label="CycE", edgecolors='none', alpha=0.5, marker="^")
plt.scatter(allthetas, resultmatrixFaure[4,:], c="#fb9a99", label="CycA", edgecolors='none', alpha=0.5, marker="*")
plt.scatter(allthetas, resultmatrixFaure[5,:], c="#a6cee3", label="p27", edgecolors='none', alpha=0.5, marker="<")
plt.scatter(allthetas, resultmatrixFaure[6,:], c="#fdbf6f", label="Cdc20", edgecolors='none', alpha=0.5, marker=">")
plt.scatter(allthetas, resultmatrixFaure[7,:], c="#ff7f00", label="Cdh1", edgecolors='none', alpha=0.5, marker="X")
plt.scatter(allthetas, resultmatrixFaure[8,:], c="#cab2d6", label="UbcH10", edgecolors='none', alpha=0.5, marker="P")
plt.scatter(allthetas, resultmatrixFaure[9,:], c="#6a3d9a", label="CycB", edgecolors='none', alpha=0.5, marker="h")

plt.legend(loc='upper right', ncol=2)
plt.grid(False)
plt.ylim((-0.05,1.05))

plt.xlabel(r'Shift in angle $\theta$ [°]', fontsize=12)
plt.ylabel("Probability of measuring the attractor state 0010100010", fontsize=12)
#plt.title(r'Deviation from unbiased single state basin size by varying angle $\theta$')
plt.title('$\\bf{b}$ Cell cycle network \n(Quantum simulation)')
ax = plt.gca()
plt.xticks(allthetas)
labelnr=0
for label in ax.xaxis.get_ticklabels():
    if labelnr % 2 != 0:
        label.set_visible(False)
    labelnr += 1

plt.axhline(y=0.5, color="black", xmin=0, xmax=180, linestyle="dashed")
plt.savefig('./FigS4/FigS4b_Faure.pdf', dpi=350)

print("Runtime = " + str(datetime.now() - starttime))

