# ColloidalMolecules
The code provided in this repository can be utilized to investigate self-assembly behavior and properties of superstructures formed by a variety of colloidal molecules. Python or Matlab codes provide the detailed description of how to generate intial configurations of various superstructures. The generated configurations can then be read by a LAMMPS input script, such as the sample script posted for the cube/octahedral system, and analyzed performing MD or BD simulations.

DNA-mediated interparticle interactions are modeled according to experimentally-validated interaction potential described in Rogers et. al, Proc. Natl. Acad. Sci. U. S. A. 2011, 108 (38), 15687–15692. Sample Tabulated interaction potentials are posted. We also note that the interaction potentials can be approximated using the mie/cut potential in LAMMPS, with gamma_rep=46, and gamma_att=23, and \epsilon representing the strength of the potential and \sigma representing a measure of the equilibrium distance between the particles. In the case of repulsive interactions, the mie/cut potential can be reproduced with a cuttof distance equal to \sigma.

A list of relevant publications is provided below:

-Parvez, N. and Zanjani, M. B. (2020) Synthetic self-limiting structures engineered with defective colloidal clusters. Adv. Funct. Mater. p. 2003317.

-Parvez, N., Rao, D. M. and Zanjani, M. B. (2019) Investigation of geometric landscape and structure–property relations for colloidal superstructures using genetic algorithm. J. Phys. Chem. B 123(34):7445–7454.

-Aryana, K., Stahley, J. B., Parvez, N., Kim, K. and Zanjani, M. B. (2019) Superstructures of multielement colloidal molecules: efficient pathways to construct reconfigurable photonic and phononic crystals. Adv. Theory Simul. p. 1800198.

-Aryana, K. and Zanjani, M. B. (2018) Diamond family of colloidal supercrystals as phononic metamaterials. J. Appl. Phys. 123(18):185103.

-Zanjani, M. B., Crocker, J. C. and Sinno, T. (2017) Self-assembly with colloidal clusters: facile crystal design using connectivity landscape analysis. Soft Matter 13(39):7098–7105.

-Zanjani, M. B., Jenkins, I. C., Crocker, J. C. and Sinno, T. (2016) Colloidal cluster assembly into ordered superstructures via engineered directional binding. ACS Nano 10(12):11280– 11289.
