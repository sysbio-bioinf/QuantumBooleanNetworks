# Quantum Boolean networks

This directory contains Python and R scripts for perfoming the analyses presented in the article *"Quantum Boolean networks for biological network analysis"* by Weidner et al.

All analysis were performed using R v4.1.0, python v3.9 and qiskit v0.36.2.
While some experiments such as those performed on the real quantum processing units are stochastic in nature, all scripts using simulators have had seeds provided.

Furthermore, additional packages are required:
For R scripts, any required packages will be directly downloaded when running the scripts if they are not yet installed.
For Python scripts, all required packages have been listed in the file `requirements.txt`.
These can be directly installed using `pipreqs` by running the following code on the command line:

```
pip install pipreqs
pip install -r requirements.txt
```

Running experiments on the IonQ quantum processing unit requires AWS account data which has to be provided by the user. Note that running such code will generate [fees](https://aws.amazon.com/de/braket/pricing/) which depends on the chosen device as well as the number of shots. More information about how to use the IonQ QPU via Amazon Braket can be found [here](https://aws.amazon.com/de/braket/).

Experiments on IBM QPUs using more than 5 qubits can only be performed by members of the [IBM Quantum Network](https://www.ibm.com/quantum/network).

