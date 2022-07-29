#Libraries: 
from qiskit import * 
from qiskit.circuit import classical_function, Int1 
 
#Regulatory functions to synthesize into circuits using Qiskit: 
@classical_function 
def g0_update(fgf8: Int1, emx2: Int1, pax6: Int1, sp8: Int1, coup_tfi: Int1) -> Int1:
	 return (fgf8 and not emx2 and sp8) 
 
@classical_function 
def g1_update(fgf8: Int1, emx2: Int1, pax6: Int1, sp8: Int1, coup_tfi: Int1) -> Int1:
	 return (not fgf8 and not pax6 and coup_tfi and not sp8) 
 
@classical_function 
def g2_update(fgf8: Int1, emx2: Int1, pax6: Int1, sp8: Int1, coup_tfi: Int1) -> Int1:
	 return (not emx2 and not coup_tfi and sp8) 
 
@classical_function 
def g3_update(fgf8: Int1, emx2: Int1, pax6: Int1, sp8: Int1, coup_tfi: Int1) -> Int1:
	 return (fgf8 and not  emx2) 
 
@classical_function 
def g4_update(fgf8: Int1, emx2: Int1, pax6: Int1, sp8: Int1, coup_tfi: Int1) -> Int1:
	 return (not fgf8 and not  sp8) 
 
