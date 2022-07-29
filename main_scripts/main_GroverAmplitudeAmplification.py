from functions import *
from qiskit import *

starttime = datetime.now()
#Perform amplitude amplification for predecessors of the 01101 attractor state in the network of Giacomantonio et al.
#Running this script required 2sec for nrshots=10000 on a MacBook Pro with a 2.3GHz Quad-Core processor and 16 GB of RAM

GiacoPath = "./28nets_rules/Giacomantonio2010.txt"
nrshots=10000

#Amplify T^inv=1 predecessors of the 01101 attractor (attractor itself + 1 extra state), using G=1 iteration
print("markedState=01101, T^inv=1, G=1:")
GiacoGroverCirc_G1_T1 = generate_groverSTGinversion_circuit(GiacoPath, nrTransitions=1, markedState=[0,1,1,0,1], G=1)
result = execute(GiacoGroverCirc_G1_T1, backend=Aer.get_backend('qasm_simulator'),
                 shots=nrshots, seed_transpiler=123, seed_simulator=123).result()
countDict_G1_T1 = result.get_counts(GiacoGroverCirc_G1_T1)
print(outputformat(countDict_G1_T1, normalizeOutput=True, sortOutput=True, digits=None))
G1_T1_solutions = countDict_G1_T1["01101"] + countDict_G1_T1["01001"]
print(str(G1_T1_solutions/nrshots) + "\n") #RESULT: 47.41% solutions for T^inv=1 and G=1

#Amplify T^inv=2 pre-predecessors of the 01101 attractor (attractor itself + 3 extra states), using G=1 iteration
print("markedState=01101, T^inv=2, G=1:")
GiacoGroverCirc_G1_T2 = generate_groverSTGinversion_circuit(GiacoPath, nrTransitions=2, markedState=[0,1,1,0,1], G=1)
result = execute(GiacoGroverCirc_G1_T2, backend=Aer.get_backend('qasm_simulator'),
                 shots=nrshots, seed_transpiler=123, seed_simulator=123).result()
countDict_G1_T2 = result.get_counts(GiacoGroverCirc_G1_T2)
print(outputformat(countDict_G1_T2, normalizeOutput=True, sortOutput=True, digits=None))
G1_T2_solutions = countDict_G1_T2["01101"] + countDict_G1_T2["01001"] + countDict_G1_T2["11001"] + countDict_G1_T2["11101"]
print(str(G1_T2_solutions/nrshots) + "\n") #RESULT: 77.68% solutions for T^inv=2 and G=1



#Further amplification is possible by increasing the number of iterations G towards G_opt, afterwards the probabilities of solution states decrease again:

#Amplify T^inv=1 predecessors of the 01101 attractor (attractor itself + 1 extra state), using G=2 iterations
print("markedState=01101, T^inv=1, G=2:")
GiacoGroverCirc_G2_T1 = generate_groverSTGinversion_circuit(GiacoPath, nrTransitions=1, markedState=[0,1,1,0,1], G=2)
result = execute(GiacoGroverCirc_G2_T1, backend=Aer.get_backend('qasm_simulator'),
                 shots=nrshots, seed_transpiler=123, seed_simulator=123).result()
countDict_G2_T1 = result.get_counts(GiacoGroverCirc_G2_T1)
print(outputformat(countDict_G2_T1, normalizeOutput=True, sortOutput=True, digits=None))
G2_T1_solutions = countDict_G2_T1["01101"] + countDict_G2_T1["01001"]
print(str(G2_T1_solutions/nrshots) + "\n") #RESULT: 47.41%->90.94% solutions for T^inv=1 and G=1->G=2 iterations

#Amplify T^inv=1 predecessors of the 01101 attractor (attractor itself + 1 extra state), using G=3 iterations
print("markedState=01101, T^inv=1, G=3:")
GiacoGroverCirc_G3_T1 = generate_groverSTGinversion_circuit(GiacoPath, nrTransitions=1, markedState=[0,1,1,0,1], G=3)
result = execute(GiacoGroverCirc_G3_T1, backend=Aer.get_backend('qasm_simulator'),
                 shots=nrshots, seed_transpiler=123, seed_simulator=123).result()
countDict_G3_T1 = result.get_counts(GiacoGroverCirc_G3_T1)
print(outputformat(countDict_G3_T1, normalizeOutput=True, sortOutput=True, digits=None))
G3_T1_solutions = countDict_G3_T1["01101"] + countDict_G3_T1["01001"]
print(str(G3_T1_solutions/nrshots) + "\n") #RESULT: 90.94%->96.12% solutions for T^inv=1 and G=2->G=3 iterations

#Amplify T^inv=1 predecessors of the 01101 attractor (attractor itself + 1 extra state), using G=4 iterations
print("markedState=01101, T^inv=1, G=4:")
GiacoGroverCirc_G4_T1 = generate_groverSTGinversion_circuit(GiacoPath, nrTransitions=1, markedState=[0,1,1,0,1], G=4)
result = execute(GiacoGroverCirc_G4_T1, backend=Aer.get_backend('qasm_simulator'),
                 shots=nrshots, seed_transpiler=123, seed_simulator=123).result()
countDict_G4_T1 = result.get_counts(GiacoGroverCirc_G4_T1)
print(outputformat(countDict_G4_T1, normalizeOutput=True, sortOutput=True, digits=None))
G4_T1_solutions = countDict_G4_T1["01101"] + countDict_G4_T1["01001"]
print(str(G4_T1_solutions/nrshots) + "\n") #RESULT: 96.12%->58.53% solutions for T^inv=1 and G=3->G=4 iterations => G=3 was optimal

#Amplify T^inv=2 pre-predecessors of the 01101 attractor (attractor itself + 3 extra states), using G=2 iterations
print("markedState=01101, T^inv=2, G=2:")
GiacoGroverCirc_G2_T2 = generate_groverSTGinversion_circuit(GiacoPath, nrTransitions=2, markedState=[0,1,1,0,1], G=2)
result = execute(GiacoGroverCirc_G2_T2, backend=Aer.get_backend('qasm_simulator'),
                 shots=nrshots, seed_transpiler=123, seed_simulator=123).result()
countDict_G2_T2 = result.get_counts(GiacoGroverCirc_G2_T2)
print(outputformat(countDict_G2_T2, normalizeOutput=True, sortOutput=True, digits=None))
G2_T2_solutions = countDict_G2_T2["01101"] + countDict_G2_T2["01001"] + countDict_G2_T2["11001"] + countDict_G2_T2["11101"]
print(str(G2_T2_solutions/nrshots) + "\n") #RESULT: 77.68%->94.60% solutions for T^inv=2 and G=1->G=2 iterations

#Amplify T^inv=2 pre-predecessors of the 01101 attractor (attractor itself + 3 extra states), using G=3 iterations
print("markedState=01101, T^inv=2, G=3:")
GiacoGroverCirc_G3_T2 = generate_groverSTGinversion_circuit(GiacoPath, nrTransitions=2, markedState=[0,1,1,0,1], G=3)
result = execute(GiacoGroverCirc_G3_T2, backend=Aer.get_backend('qasm_simulator'),
                 shots=nrshots, seed_transpiler=123, seed_simulator=123).result()
countDict_G3_T2 = result.get_counts(GiacoGroverCirc_G3_T2)
print(outputformat(countDict_G3_T2, normalizeOutput=True, sortOutput=True, digits=None))
G3_T2_solutions = countDict_G3_T2["01101"] + countDict_G3_T2["01001"] + countDict_G3_T2["11001"] + countDict_G3_T2["11101"]
print(str(G3_T2_solutions/nrshots) + "\n") #RESULT: 94.60%->33.62% solutions for T^inv=2 and G=2->G=3 iterations => G=2 was optimal

#Conclusion: G_opt=3 for M/N=2/32 solutions and G_opt=2 for M/N=4/32 solutions, as given by equation for G_opt


print("Try amplification of (nonexistent) predecessors of a GoE state, such as 01100:")
#Amplify T^inv=1 predecessors of the 01100 GoE state, i.e. M=0, using G=1 iterations
print("markedState=01100 (GoE state), T^inv=1, G=1:")
GiacoGroverCirc_G1_T1_GoE = generate_groverSTGinversion_circuit(GiacoPath, nrTransitions=1, markedState=[0,1,1,0,0], G=1)
result = execute(GiacoGroverCirc_G1_T1_GoE, backend=Aer.get_backend('qasm_simulator'),
                 shots=nrshots, seed_transpiler=123, seed_simulator=123).result()
countDict_G1_T1_GoE = result.get_counts(GiacoGroverCirc_G1_T1_GoE)
print(outputformat(countDict_G1_T1_GoE, normalizeOutput=True, sortOutput=True, digits=None))
#RESULT: Same uniform superposition returned as output that was given as input, no amplification has taken place

print("Runtime = " + str(datetime.now() - starttime))
