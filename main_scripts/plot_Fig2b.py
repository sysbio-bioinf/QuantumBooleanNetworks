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
#Running this script required 2min47sec for shotnr=10000 on a MacBook Pro with a 2.3GHz Quad-Core processor and 16 GB of RAM

GiacoPath = "./28nets_rules/Giacomantonio2010.txt"

# OUTER PLOT:
labelints = list(range(32))
labelbitstrings = ['{0:05b}'.format(x) for x in labelints] #5 genes in network, 5 bits shown for state labels

#RUN STATE TRANSITION CIRCUITS
shotnr = 10000 #Number of measurements performed for each T
#T=0
GiacoExactTransition0 = allparam_exact_multiTransition_synchronous(rulestxt=GiacoPath, Tmax=0, nrshots=shotnr, normalizeOutput=True, seed_transpiler=123, seed_simulator=123)
print(GiacoExactTransition0)
#T=1
GiacoExactTransition1 = allparam_exact_multiTransition_synchronous(rulestxt=GiacoPath, Tmax=1, nrshots=shotnr, normalizeOutput=True, seed_transpiler=123, seed_simulator=123)
print(GiacoExactTransition1)
#T=2
GiacoExactTransition2 = allparam_exact_multiTransition_synchronous(rulestxt=GiacoPath, Tmax=2, nrshots=shotnr, normalizeOutput=True, seed_transpiler=123, seed_simulator=123)
print(GiacoExactTransition2)
#T=3
GiacoExactTransition3 = allparam_exact_multiTransition_synchronous(rulestxt=GiacoPath, Tmax=3, nrshots=shotnr, normalizeOutput=True, seed_transpiler=123, seed_simulator=123)
print(GiacoExactTransition3)
#T=4
GiacoExactTransition4 = allparam_exact_multiTransition_synchronous(rulestxt=GiacoPath, Tmax=4, nrshots=shotnr, normalizeOutput=True, seed_transpiler=123, seed_simulator=123)
print("Final attractor distribution")
print(GiacoExactTransition4)

if pickledump:
    pickle.dump(GiacoExactTransition0, open("./Fig2b/T0dictionary.pkl", "wb"))
    pickle.dump(GiacoExactTransition1, open("./Fig2b/T1dictionary.pkl", "wb"))
    pickle.dump(GiacoExactTransition2, open("./Fig2b/T2dictionary.pkl", "wb"))
    pickle.dump(GiacoExactTransition3, open("./Fig2b/T3dictionary.pkl", "wb"))
    pickle.dump(GiacoExactTransition4, open("./Fig2b/T4dictionary.pkl", "wb"))
else:
    GiacoExactTransition0 = pickle.load(open("./Fig2b/T0dictionary.pkl", "rb"))
    GiacoExactTransition1 = pickle.load(open("./Fig2b/T1dictionary.pkl", "rb"))
    GiacoExactTransition2 = pickle.load(open("./Fig2b/T2dictionary.pkl", "rb"))
    GiacoExactTransition3 = pickle.load(open("./Fig2b/T3dictionary.pkl", "rb"))
    GiacoExactTransition4 = pickle.load(open("./Fig2b/T4dictionary.pkl", "rb"))


# INSET PLOT: Load theta variation data
print("Begin calculation for initial theta variation in Giacomantonio network... ")
FullSyncCirc = synthesizeFullNetworkUpdateCircuit(rulestxt=GiacoPath, update="synchronous")
N = FullSyncCirc.num_qubits//2
Tmax_Giaco = 4
allthetas = list(range(0,190, 10)) #0,...,180
thetainitvec = list(np.repeat(90, N))

resultmatrix = np.zeros((N, len(allthetas)), dtype=float)

#INSET: Calculation of attractor probability shift by changing angle theta:

for g in range(N):
    print("Gene " + str(g))
    for thetanr in range(len(allthetas)):
        thetashift = allthetas[thetanr]
        print("thetashift " + str(thetashift))
        # Change to only vary thetainitvec in gene g
        thetainitvec_g = copy.deepcopy(thetainitvec)
        thetainitvec_g[g] = thetashift

        Transition = allparam_exact_multiTransition_synchronous(rulestxt=GiacoPath, Tmax=Tmax_Giaco,
                                                         nrshots=shotnr, initActivities=thetainitvec_g,
                                                         thetatype="angle", normalizeOutput=True,
                                                         seed_transpiler=123, seed_simulator=123)
        #was calling allparam_multiTransition_synchronous before -> changed to _exact_, need the other one only for Faure circuits

        if "10010" in Transition.keys():
            resultmatrix[g,thetanr] = Transition["10010"]
        else:
            print("Attractor not measured!")
            resultmatrix[g,thetanr] = 0

        #Save results
        if pickledump:
            np.save(file="./Fig2b/GiacomantonioShiftByTheta.npy", arr=resultmatrix)


print("Begin plotting Giacomantonio network... ")

#Load previous results
if pickledump == False:
    resultmatrix = np.load(file="./Fig2b/GiacomantonioShiftByTheta.npy")


# PLOT GENERATION: barplot with convergence to attractor states in quantum circuit, shift by varying initial state with R_y(theta) gates as inset
# set width of bars
barWidth = 0.15
# set heights of bars (y values to plot)
bars0 = probdict2list(GiacoExactTransition0)
bars1 = probdict2list(GiacoExactTransition1)
bars2 = probdict2list(GiacoExactTransition2)
bars3 = probdict2list(GiacoExactTransition3)
bars4 = probdict2list(GiacoExactTransition4)
# Set position of bar on X axis
r0 = np.arange(len(bars0))
r1 = [x + barWidth for x in r0]
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]
r4 = [x + barWidth for x in r3]

fig, ax = plt.subplots()
fig.set_size_inches(w=10, h=5.5)
ax.bar(r0, bars0, color='#000000', width=barWidth, edgecolor='white', label='T=0')
ax.bar(r1, bars1, color='#557f2d', width=barWidth, edgecolor='white', label='T=1')
ax.bar(r2, bars2, color='#002eff', width=barWidth, edgecolor='white', label='T=2')
ax.bar(r3, bars3, color='#800080', width=barWidth, edgecolor='white', label='T=3')
ax.bar(r4, bars4, color='#ff0000', width=barWidth, edgecolor='white', label='T=4')
ax.set_ylabel("Probability", fontsize=18)
# Add xticks on the middle of the group bars
ax.set_xlabel('State', fontsize=18)
ax.set_title("$\\bf{b}$ Convergence to attractor states in a multi-transition circuit", fontsize=15)
ax.set_yticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
ax.legend(framealpha=1, loc="upper left")
ax.grid(color='gray', linestyle='dashed', axis='y', alpha=0.5)

ax.set_xticks(list(range(0,32)))
ax.xaxis.set_ticklabels(["00000", "00001", "00010", "00011", "00100",
               "00101", "00110", "00111", "01000", "01001",
               "01010", "01011", "01100", "$\\bf{01101}$", "01110",
               "01111", "10000", "10001", "$\\bf{10010}$", "10011",
               "10100", "10101", "10110", "10111", "11000",
               "11001", "11010", "11011", "11100", "11101",
               "11110", "11111"], rotation=75, ha="center")

axIns = ax.inset_axes([0.65, 0.5, 0.3, 0.3])
axIns.scatter(allthetas, resultmatrix[0,:], c="#e31a1c", label="Fgf8", edgecolors='none', alpha=0.5, marker="o")
axIns.scatter(allthetas, resultmatrix[1,:], c="#cab2d6", label="Emx2", edgecolors='none', alpha=0.75, marker="s")
axIns.scatter(allthetas, resultmatrix[2,:], c="#33a02c", label="Pax6", edgecolors='none', alpha=0.5, marker="v")
axIns.scatter(allthetas, resultmatrix[3,:], c="#1f78b4", label="Sp8", edgecolors='none', alpha=0.75, marker="^")
axIns.scatter(allthetas, resultmatrix[4,:], c="#fdbf6f", label="Coup-tfi", edgecolors='none', alpha=1, marker="*")

axIns.legend(loc='lower right', fontsize=5, ncol=2)
axIns.grid(False)
axIns.set_title('Variation in probability of attractor $\\bf{10010}$\ndepending on initial state bias', fontsize=9)
axIns.set_ylim((0.55,1.05))
axIns.set_ylabel("Probability", fontsize=9)
axIns.set_xlabel(r'Shift in angle $\theta$ [°]', fontsize=9)
axIns.set_xticks([0,45,90,135,180])
axIns.set_yticks([0.6,0.7,0.8,0.9,1.0])
axIns.axhline(y=0.875, color="red", xmin=0, xmax=180, linestyle="dashed")

plt.savefig('./Fig2b/Fig2b.pdf', dpi=350)

print("Runtime = " + str(datetime.now() - starttime))

