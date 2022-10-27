# G-FE-G-tunneling

A code that calculates the junction electrostatics and I-V curve for ferroelectrically modulated tunnel junctions, uploaded here as a supplement to the article and master thesis project I did with prof. Eran Sela.

The main .m file is "numeric_qc_debug", in which the parameter models of the junction should be specified. 
Based on the that, the code calculates the junction energetics by minimizing the energy functional for a given gate and bias voltage. After the energetics is found, the Bardeen transfer Hamiltonian approach is used to calculate the I-V characteristic of the system.

More details could be found in the relavant article [here](https://doi.org/10.48550/arXiv.2206.13249) or [here](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.106.144110).

![Device_new_figure](https://user-images.githubusercontent.com/109433383/179699030-6a8537bb-0018-41b1-bb3a-c67662c50191.svg)
