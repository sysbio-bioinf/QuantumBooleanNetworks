# Quantum Boolean networks

This directory contains Python and R scripts for perfoming the analyses presented in the article *"Quantum Boolean networks for molecular network analysis"* by Weidner et al.

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

The required R packages are:
- gtools v3.9.2
- BoolNet v2.1.5
- poweRlaw v0.70.6
- ggplot2 v3.3.5
- scales 1.1.1
These packages will be downloaded when running the respective R scripts if they are not yet installed.

Running experiments on the IonQ quantum processing unit requires AWS account data which has to be provided by the user. Note that running such code will generate [fees](https://aws.amazon.com/de/braket/pricing/) which depends on the chosen device as well as the number of shots. More information about how to use the IonQ QPU via Amazon Braket can be found [here](https://aws.amazon.com/de/braket/).

Experiments on IBM QPUs using more than 5 qubits can only be performed by members of the [IBM Quantum Network](https://www.ibm.com/quantum/network).

In the following the run times for performing all analyses and generating the figures shown in the article are listed for a MacBook Pro with a 2.3GHz Quad-Core processor 16 GB of RAM.

<u>*Run times for analyses in the main manuscript:*</u>
- plot_Fig2b.py: 2min 47sec
- main_GroverAmplitudeAmplification.py: 2sec
- plot_Fig3.py: 2min 31sec
- plot_Fig4.py: This analysis depends on results from cloud-accessible quantum processing units. The time to obtain these results depends on the availability of the devices and their queue at the time of task submission. Once this data has been received, the script has a run time of 0.6sec.

<u>*Run times for analyses in the supplementary material:*</u>
- scalefree_powerlaw_checking.R: 1.06min
- plot_FigS1.py: 1.3sec
- plot_FigS4_classical.R: 1.49h
- plot_FigS4_quantum.py: 10min 18sec
- TableS2.R: 0.02sec
- TableS2.py: 0.6sec
- FigS7_generate_tt_csv.R: 7.3h
- plot_FigS7: 14.4min if transition tables from FigS7_generate_tt_csv.R are available.
- plot_FigS8: This analysis depends on results from cloud-accessible quantum processing units. The time to obtain these results depends on the availability of the devices and their queue at the time of task submission. Once this data has been received, the script has a run time of 0.6sec.

<u>*Installation time:*</u>
On the same MacBook Pro, the setup of a new virtual environment followed by the installation of the packages specified in the `requirements.txt` file required 38.5 sec.

<u>*Demo:*</u>
The main_scripts folder contains a demo script for performing dynamic analyses of a small QBN for the mammalian cortical area development network, such as state transitions, amplitude amplification and quantum counting.
