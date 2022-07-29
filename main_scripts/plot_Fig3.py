from functions import *
import numpy as np
import matplotlib.pyplot as plt
import pickle
from itertools import starmap
from operator import mul

starttime = datetime.now()

GiacoPath = "./28nets_rules/Giacomantonio2010.txt"
#Averages of Distributions to be plotted as dashed vertical lines
def avgOfDistr(distr):
    #State is already normalised, need to sum up all key*val, key is nrResults of counting, val is prob.weight
    avg = sum(starmap(mul, distr.items())) / sum(distr.values())
    return(avg)

#Note that measurement is stochastic, results will differ slightly on each run
#Load dictionaries from pkl files to reproduce the plot as shown in the article
#Running this script required 2min31sec for nrshots=1000 on a MacBook Pro with a 2.3GHz Quad-Core processor and 16 GB of RAM

#RUN QUANTUM COUNTING ALGORITHM:
pickledump = False #Should results of this run be saved (True), else load previously stored results (False)
#Run algorithm for T^inv=1 inverse transition from 01101 attractor state with various register lengths
r4_01101_Tinv1 = QuantumCountingAlgo(rulestxt=GiacoPath, nrTransitions=1, markedState=[0,1,1,0,1], r_registerLen=4, nrshots = 1000,
                                     seed_transpiler=123, seed_simulator=123)
print("Results of 1000 measurements from counting T^inv=1 predecessors of 01101 attractor state with a register of r=4 readout qubits:")
print(r4_01101_Tinv1)
avg_r4_01101_Tinv1 = avgOfDistr(r4_01101_Tinv1)
print(avg_r4_01101_Tinv1)
print("\n")

r5_01101_Tinv1 = QuantumCountingAlgo(rulestxt=GiacoPath, nrTransitions=1, markedState=[0,1,1,0,1], r_registerLen=5, nrshots = 1000,
                                     seed_transpiler=123, seed_simulator=123)
print("Results of 1000 measurements from counting T^inv=1 predecessors of 01101 attractor state with a register of r=5 readout qubits:")
print(r5_01101_Tinv1)
avg_r5_01101_Tinv1 = avgOfDistr(r5_01101_Tinv1)
print(avg_r5_01101_Tinv1)
print("\n")

r6_01101_Tinv1 = QuantumCountingAlgo(rulestxt=GiacoPath, nrTransitions=1, markedState=[0,1,1,0,1], r_registerLen=6, nrshots = 1000,
                                     seed_transpiler=123, seed_simulator=123)
print("Results of 1000 measurements from counting T^inv=1 predecessors of 01101 attractor state with a register of r=6 readout qubits:")
print(r6_01101_Tinv1)
avg_r6_01101_Tinv1 = avgOfDistr(r6_01101_Tinv1)
print(avg_r6_01101_Tinv1)
print("\n")

#Run algorithm for T^inv=2 inverse transition from 01101 attractor state with various register lengths
r4_01101_Tinv2 = QuantumCountingAlgo(rulestxt=GiacoPath, nrTransitions=2, markedState=[0,1,1,0,1], r_registerLen=4, nrshots = 1000,
                                     seed_transpiler=123, seed_simulator=123)
print("Results of 1000 measurements from counting T^inv=2 predecessors of 01101 attractor state with a register of r=4 readout qubits:")
print(r4_01101_Tinv2)
avg_r4_01101_Tinv2 = avgOfDistr(r4_01101_Tinv2)
print(avg_r4_01101_Tinv2)
print("\n")

r5_01101_Tinv2 = QuantumCountingAlgo(rulestxt=GiacoPath, nrTransitions=2, markedState=[0,1,1,0,1], r_registerLen=5, nrshots = 1000,
                                     seed_transpiler=123, seed_simulator=123)
print("Results of 1000 measurements from counting T^inv=2 predecessors of 01101 attractor state with a register of r=5 readout qubits:")
print(r5_01101_Tinv2)
avg_r5_01101_Tinv2 = avgOfDistr(r5_01101_Tinv2)
print(avg_r5_01101_Tinv2)
print("\n")

r6_01101_Tinv2 = QuantumCountingAlgo(rulestxt=GiacoPath, nrTransitions=2, markedState=[0,1,1,0,1], r_registerLen=6, nrshots = 1000,
                                     seed_transpiler=123, seed_simulator=123)
print("Results of 1000 measurements from counting T^inv=2 predecessors of 01101 attractor state with a register of r=6 readout qubits:")
print(r6_01101_Tinv2)
avg_r6_01101_Tinv2 = avgOfDistr(r6_01101_Tinv2)
print(avg_r6_01101_Tinv2)
print("\n")

if pickledump:
    pickle.dump(r4_01101_Tinv1, open("./Fig3/r4_01101_Tinv1.pkl", "wb"))
    pickle.dump(r5_01101_Tinv1, open("./Fig3/r5_01101_Tinv1.pkl", "wb"))
    pickle.dump(r6_01101_Tinv1, open("./Fig3/r6_01101_Tinv1.pkl", "wb"))
    pickle.dump(r4_01101_Tinv2, open("./Fig3/r4_01101_Tinv2.pkl", "wb"))
    pickle.dump(r5_01101_Tinv2, open("./Fig3/r5_01101_Tinv2.pkl", "wb"))
    pickle.dump(r6_01101_Tinv2, open("./Fig3/r6_01101_Tinv2.pkl", "wb"))
else:
    #ALTERNATIVELY, LOAD DATA OF PREVIOUS RUN USED FOR PLOTTING:
    r4_01101_Tinv1 = pickle.load(open("./Fig3/r4_01101_Tinv1.pkl", "rb")) #Avg=2.666
    r5_01101_Tinv1 = pickle.load(open("./Fig3/r5_01101_Tinv1.pkl", "rb")) #Avg=3.483
    r6_01101_Tinv1 = pickle.load(open("./Fig3/r6_01101_Tinv1.pkl", "rb")) #Avg=2.276
    r4_01101_Tinv2 = pickle.load(open("./Fig3/r4_01101_Tinv2.pkl", "rb")) #Avg=5.225
    r5_01101_Tinv2 = pickle.load(open("./Fig3/r5_01101_Tinv2.pkl", "rb")) #Avg=5.191
    r6_01101_Tinv2 = pickle.load(open("./Fig3/r6_01101_Tinv2.pkl", "rb")) #Avg=4.693
    avg_r4_01101_Tinv1 = avgOfDistr(r4_01101_Tinv1)
    avg_r5_01101_Tinv1 = avgOfDistr(r5_01101_Tinv1)
    avg_r6_01101_Tinv1 = avgOfDistr(r6_01101_Tinv1)
    avg_r4_01101_Tinv2 = avgOfDistr(r4_01101_Tinv2)
    avg_r5_01101_Tinv2 = avgOfDistr(r5_01101_Tinv2)
    avg_r6_01101_Tinv2 = avgOfDistr(r6_01101_Tinv2)


#GENERATE BAR PLOT
fig = plt.figure()
ax = fig.add_subplot()

xticklabels = list(range(33))
xticklabels[2] = "$\\bf{" + "2" + "}$"
xticklabels[4] = "$\\bf{" + "4" + "}$"
plt.xticks(list(range(33)), xticklabels) #0-32
plt.yticks([x/10 for x in list(range(0,11))])
plt.ylim((0,1))
plt.xlabel("Number of predecessor states M", fontsize=18)
plt.ylabel("Probability", fontsize=18)
plt.title("Counting predecessors of 01101 attractor state using different \n numbers of transitions T$^{inv}$ and readout register sizes $r$", fontsize=15)
ax.yaxis.grid()

plt.bar(r4_01101_Tinv1.keys(), r4_01101_Tinv1.values(), width=1, color='gold', alpha=0.5, label="r=4 T$^{inv}$ = 1")
plt.bar(r5_01101_Tinv1.keys(), r5_01101_Tinv1.values(), width=1, color='darkorange', alpha=0.5, label="r=5, T$^{inv}$ = 1")
plt.bar(r6_01101_Tinv1.keys(), r6_01101_Tinv1.values(), width=1, color='maroon', alpha=0.5, label="r=6, T$^{inv}$ = 1")
plt.bar(r4_01101_Tinv2.keys(), r4_01101_Tinv2.values(), width=1, color='lightblue', alpha=0.5, label="r=4, T$^{inv}$ = 2")
plt.bar(r5_01101_Tinv2.keys(), r5_01101_Tinv2.values(), width=1, color='blue', alpha=0.5, label="r=5, T$^{inv}$ = 2")
plt.bar(r6_01101_Tinv2.keys(), r6_01101_Tinv2.values(), width=1, color='darkblue', alpha=0.5, label="r=6, T$^{inv}$ = 2")

plt.legend(framealpha=1)
plt.tight_layout()

plt.axvline(avg_r4_01101_Tinv1, color='gold', linestyle='dashed', linewidth=2, alpha=1)
plt.axvline(avg_r5_01101_Tinv1, color='darkorange', linestyle='dashed', linewidth=2, alpha=1)
plt.axvline(avg_r6_01101_Tinv1, color='maroon', linestyle='dashed', linewidth=2, alpha=1)
plt.axvline(avg_r4_01101_Tinv2, color='lightblue', linestyle='dashed', linewidth=2, alpha=1)
plt.axvline(avg_r5_01101_Tinv2, color='blue', linestyle='dashed', linewidth=2, alpha=1)
plt.axvline(avg_r6_01101_Tinv2, color='darkblue', linestyle='dashed', linewidth=2, alpha=1)

#plt.show()
fig.set_size_inches(10, 5)
plt.savefig('./Fig3/Fig3.pdf', dpi=350)
plt.clf()

print("Runtime = " + str(datetime.now() - starttime))

#Test: Show influence of GoE states:
print("Predecessors of 00100, show influence of GoE states:\n")
r6_onestep00100 = QuantumCountingAlgo(rulestxt=GiacoPath, nrTransitions=1, markedState=[0,0,1,0,0], r_registerLen=6, nrshots = 1000,
                                      seed_transpiler=123, seed_simulator=123)
print("Influence of GoE states on predecessors of 00100 = 2 predecessors, 1 of these is a GoE state:")
print(outputformat(r6_onestep00100, digits=3))

r6_twostep00100 = QuantumCountingAlgo(rulestxt=GiacoPath, nrTransitions=2, markedState=[0,0,1,0,0], r_registerLen=6, nrshots = 1000,
                                      seed_transpiler=123, seed_simulator=123)
print("Influence of GoE states on pre-predecessors of 00100 = 4 pre-predecessors:")
print(outputformat(r6_twostep00100, digits=3))
