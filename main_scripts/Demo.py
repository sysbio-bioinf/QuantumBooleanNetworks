#Import necessary functions, set number of measurements to perform for mammalian cortical area development network
from functions import *
GiacoPath = "./28nets_rules/Giacomantonio2010.txt"
nrshots=10000

#Print a circuit for performing a synchronous state transition:
syncTransitionCircuit = synthesizeFullNetworkUpdateCircuit(GiacoPath)
print(syncTransitionCircuit)

#Perform T=4 forward state transitions in the network.
#Perform shotnr repeated measurements to obtain the attractor distribution.
GiacoExactTransition4 = allparam_exact_multiTransition_synchronous(rulestxt=GiacoPath, Tmax=4, nrshots=nrshots, normalizeOutput=True, seed_transpiler=123, seed_simulator=123)
print("Final attractor distribution:")
print(GiacoExactTransition4)

#Perform amplitude amplification on the chosen marked state, do T^inv inverted transitions.
#Choose the number G of iterations of the Grover operator to apply:
print("markedState=01101, T^inv=1, G=3:")
GiacoGroverCirc_G3_T1 = generate_groverSTGinversion_circuit(GiacoPath, nrTransitions=1, markedState=[0,1,1,0,1], G=3)
result = execute(GiacoGroverCirc_G3_T1, backend=Aer.get_backend('qasm_simulator'),
                 shots=nrshots, seed_transpiler=123, seed_simulator=123).result()
countDict_G3_T1 = result.get_counts(GiacoGroverCirc_G3_T1)
print(outputformat(countDict_G3_T1, normalizeOutput=True, sortOutput=True, digits=None))
G3_T1_solutions = countDict_G3_T1["01101"] + countDict_G3_T1["01001"]
print(str(G3_T1_solutions/nrshots) + "\n") #RESULT: 96.12% solutions for T^inv=1 and G_opt = 3 iterations


print("markedState=01101, T^inv=2, G=2:")
GiacoGroverCirc_G2_T2 = generate_groverSTGinversion_circuit(GiacoPath, nrTransitions=2, markedState=[0,1,1,0,1], G=2)
result = execute(GiacoGroverCirc_G2_T2, backend=Aer.get_backend('qasm_simulator'),
                 shots=nrshots, seed_transpiler=123, seed_simulator=123).result()
countDict_G2_T2 = result.get_counts(GiacoGroverCirc_G2_T2)
print(outputformat(countDict_G2_T2, normalizeOutput=True, sortOutput=True, digits=None))
G2_T2_solutions = countDict_G2_T2["01101"] + countDict_G2_T2["01001"] + countDict_G2_T2["11001"] + countDict_G2_T2["11101"]
print(str(G2_T2_solutions/nrshots) + "\n") #RESULT: 94.60% solutions for T^inv=2 and G_opt = 2 iterations

#Perform the quantum counting algorithm to estimate the number of solutions M of all T^inv=nrTransitions predecessors
#of the marked state. Accuracy of the results depends on the size of the readout register r_registerLen
nrshots = 1000
r6_01101_Tinv2 = QuantumCountingAlgo(rulestxt=GiacoPath, nrTransitions=2, markedState=[0,1,1,0,1], r_registerLen=6, nrshots = nrshots,
                                     seed_transpiler=123, seed_simulator=123)
print("Results of 1000 measurements from counting T^inv=2 predecessors of 01101 attractor state with a register of r=6 readout qubits:")
print(r6_01101_Tinv2)
avg_r6_01101_Tinv2 = avgOfDistr(r6_01101_Tinv2)
print(avg_r6_01101_Tinv2)
#Given the STG of the classical BN, the basin size of the 01101 attractor is 4/32.

