# SLAM-3D-Extension

01 Aug 2022 - FDM check done for lattice vector derivatives ✌️

04 Sep 2022 - ShellModel / PointCharge electrostatic calculations : energy, 1st derivatives, strain derivatives, lattice derivatives done.

29 Sep 2022 - Integrating SLAM LP atoms, Real space LP - LP (i.e., density - density) interaction recipes 
              (FDM check done - energy / 1st derivatives)

02 Nov 2022 - All the required integrals for the LonePair related interactions were coded (i.e., their mathematical formula) and confirmed by using
              the finite-difference-method.
              
  * Real space integrals       - real: energy / 1st / 2nd derivatives - essential to compute LP-LP interactions + atomic forces
  * Reciprocal space integrals - imag: energy / 1st derivatives
              
  This means, all the recipes has been coded, the left tasks are assembling them up to get energy/forces/strain derivatives of the lattice
  
  * Real space part has been updated so far, and the self-consistent field method was applied to see if the system shows convergence.
    (Tested for various configurations, more than one LonePair plus shell and rigid-ion in the unitcell)
  
  * some acceleration is required though, the real space integrals take ~ 600(s) to compute, therefore careful adjustments of the Ewald
    summation parameters (known as sigma) is required and also some technical support, e.g., MPI/OpenMP, implementation will be added in the next update.
