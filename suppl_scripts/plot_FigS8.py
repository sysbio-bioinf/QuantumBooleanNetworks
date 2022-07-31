from functions import *
import pickle
import matplotlib.pyplot as plt

starttime = datetime.now()
#Running this script required 0.6sec on a MacBook Pro with a 2.3GHz Quad-Core processor and 16 GB of RAM

#Plot Figure S8 comparing the results of a single state transition in the network of Giacomantonio et al.
#Compare a noiseless state vector simulator to results from IBMQ Toronto (both mitigated and unmitigated) as well as a noisy mock backend


GiacoPath = "./28nets_rules/Giacomantonio2010.txt"

#Distribution of 1/32 probabilities across all 2^5 = 32 states
uniformDistribution = pickle.load(open("./FigS8/uniformDistribution_5qubits.pkl", "rb"))

#Simulation for ideal results, without transpilation, without noise
#Noiseless_1kshots = allparam_exact_multiTransition_synchronous(GiacoPath, Tmax=1, nrshots=1000, seed_transpiler=123, seed_simulator=123)
#pickle.dump(Noiseless_1kshots, open("./FigS8/Noiseless_1kshots.pkl", "wb"))
Noiseless_1kshots = pickle.load(open("./FigS8/Noiseless_1kshots.pkl", "rb"))

#Real IBMQ Toronto QPU + dynamical decoupling + readout error mitigation:
Toronto_dd_mitigated = pickle.load(open("./FigS8/Toronto_dd_mitigated.pkl", "rb"))

#Real IBMQ Toronto, unmitigated:
Toronto_nodd_nomitigation = pickle.load(open("./FigS8/Toronto_nodd_nomitigation.pkl", "rb"))

#Fake IBMQ Toronto QPU - mock backend, using same transpilation and simulator seeds as noiseless simulator:
from qiskit.test.mock import FakeToronto
TorontoBackend = FakeToronto()
#MockBackend_NoisyFakeToronto = allparam_exact_multiTransition_synchronous(GiacoPath, Tmax=1, nrshots=1000,
#                                                                          seed_transpiler=123, seed_simulator=123,
#                                                                          addNoise=True, backend=TorontoBackend, transpileCircuit=True)
MockBackend_NoisyFakeToronto = pickle.load(open("./FigS8/MockBackend_NoisyFakeToronto.pkl", "rb"))

#Get normalized Hellinger fidelities
#Ideal = Noiseless simulator
normedFidelity = normalized_fidelity_Lubinski(idealDistr=Noiseless_1kshots, outputDistr=Noiseless_1kshots, uniformDistr=uniformDistribution)
print("Normalised Hellinger fidelity, ideal=Noiseless simulator, output=Noiseless simulator = " + str(round(normedFidelity,3))) #1.0
normedFidelity = normalized_fidelity_Lubinski(idealDistr=Noiseless_1kshots, outputDistr=Toronto_dd_mitigated, uniformDistr=uniformDistribution)
print("Normalised Hellinger fidelity, ideal=Noiseless simulator, output=IBMQ Toronto (mitigated) = " + str(round(normedFidelity,3))) #0.196
normedFidelity = normalized_fidelity_Lubinski(idealDistr=Noiseless_1kshots, outputDistr=Toronto_nodd_nomitigation, uniformDistr=uniformDistribution)
print("Normalised Hellinger fidelity, ideal=Noiseless simulator, output=IBMQ Toronto (unmitigated) = " + str(round(normedFidelity,3))) #0.114
normedFidelity = normalized_fidelity_Lubinski(idealDistr=Noiseless_1kshots, outputDistr=MockBackend_NoisyFakeToronto, uniformDistr=uniformDistribution)
print("Normalised Hellinger fidelity, ideal=Noiseless simulator, output=FakeToronto = " + str(round(normedFidelity,3))) #0.250
print("\n")

# Generate barplot

# set width of bars
barWidth = 0.25
# set heights of bars (y values to plot)
bars0 = probdict2list(Noiseless_1kshots)
bars1 = probdict2list(MockBackend_NoisyFakeToronto)
bars2 = probdict2list(Toronto_dd_mitigated)
bars3 = probdict2list(Toronto_nodd_nomitigation)
r0 = np.arange(len(bars0))
r1 = [x + barWidth for x in r0]
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]

# Make the plot
labelints = list(range(32))
labelbitstrings = ['{0:05b}'.format(x) for x in labelints] #5 bits shown in total
figure = plt.gcf()
plt.bar(r0, bars0, color='#000000', width=barWidth, edgecolor='white', label='Noiseless simulator')
plt.bar(r1, bars1, color='#0fff37', width=barWidth, edgecolor='white', label='FakeToronto')
plt.bar(r2, bars2, color='#003dca', width=barWidth, edgecolor='white', label='Error mitigated QPU')
plt.bar(r3, bars3, color='#0fedff', width=barWidth, edgecolor='white', label='Unmitigated QPU')

plt.ylabel("Probability", fontsize=18)
# Add xticks on the middle of the group bars
plt.xlabel('State', fontsize=18)

remainingstateindices = [0,4,8,9,13,16,18] #Indices of states which have non-zero weight in exact transition, highlight in bold font
for a in remainingstateindices:
    labelbitstrings[a] = "$\\bf{" + labelbitstrings[a] + "}$"

plt.xticks([r + barWidth for r in range(len(bars0))], labelbitstrings)
plt.xticks(rotation=75, ha="center")
plt.title("IBMQ Toronto error comparisons", fontsize=15)

plt.grid(color='gray', linestyle='dashed', axis='y')
plt.yticks([0.1,0.2,0.3,0.4,0.5])

# Create legend & Show graphic
plt.legend(framealpha=1)
#plt.show()
figure.set_size_inches(w=10, h=5)
plot_margin = 0.25
plt.tight_layout()
plt.savefig('./FigS8/FigS8.pdf', dpi=350)

print("Runtime = " + str(datetime.now() - starttime))
