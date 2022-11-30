# SLAM-3D-Extension - Development Log

01 Aug 2022 - FDM check done for lattice vector derivatives.

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
  
  * Some acceleration is required though, the real space integrals take ~ 600(s) to compute, therefore careful adjustments of the Ewald
    summation parameters (known as sigma) is necessary and some technical supports, e.g., MPI/OpenMP, implementation will be added in the next update.

30 Nov 2022 - Energy Calculation

  1. Lattice energy (Coulomb or electrostatic) calculation is now available.
  - Various cell energy calculations have been done, from cubic rocksalt-like to monoclinic cell. Calculated energy values were compared with GULP and verified
    the model energy is reliable. List of tested cells: alpha-PbO, perovskite ABX3, MgO-like cell/supercell (2x2x2) - tilting/stretching of lattice vectors / dislocation of atoms were also tested.
    
  2. More than one lone pair atoms (or cations) within a cell, current version of code carries out the Self-Consistent Field (SCF) calculations.
  - There is no acceleration method applied, however, simple damping or DIIS algorithm may be included in the next update.
    
  3. Performance Update - The efficiency of the code has been remarkably improved.
  - Benchmark: a single SCF calculation wtime, ~20s -> ~0.7s (alpha-PbO cell tested).
  - Balancing in the Ewald summation, real/reciprocal, was done (current sigma parameter : 1.123 Ansgtrom), without loosing the required accuracy, 10E-8 (eV).
  - OpenMP was applied to the raidal integrals (tested on UCL-RC, MMM Young, single node, 36 thread).
  - Radial integration method has been changed to Gauss-Quad from Trapezoidal.
  - Memoization done, removing out all the unnecessary repeated calcultions.
