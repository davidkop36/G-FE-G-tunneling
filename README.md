# G-FE-G-tunneling

A code that calculates the junction electrostatics and I-V curve for ferroelectrically modulated tunnel junctions. 

The main .m file is "numeric_qc_debug", in which the parameter models of the junction should be specified. 
Based on the that, the code calculates the junction energetics by minimizing the energy functional for a given gate and bias voltage. After the energetics is found, the Bardeen transfer Hamiltonian approach is used to calculate the I-V characteristic of the system.


![Device_new_figure](https://user-images.githubusercontent.com/109433383/179699030-6a8537bb-0018-41b1-bb3a-c67662c50191.svg)
