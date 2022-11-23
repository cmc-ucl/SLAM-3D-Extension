#include <iostream>
//#include <mpi.h>

#include "Cell.hpp"

#define TARGET_INPUT "geo.txt"

int main()
{
	

	Cell c(TARGET_INPUT);
	c.ShowBasicCellInfo();

	// Calculate Core / Core-Shell Contribution
	c.CalcCoulombEnergy();
	c.CalcCoulombDerivative();

	// Calculate LonePair Contribution
	c.CalcLonePairCoulombEnergy();			// This is not globally working when there are now LP cations ... need modification !!! 27 Oct 2022
	c.CalcLonePairCoulombDerivative();


	// This has to be called after StrainDerivatives are ready
	c.CalcLatticeDerivative();
	c.ShowEnergyDerivative();

	c.Finalise();
	return 0;
}
