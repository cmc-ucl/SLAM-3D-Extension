#include <iostream>
#include "cell.hpp"

#define TARGET_INPUT "geo.txt"

int main()
{
	Cell c(TARGET_INPUT);
	c.ShowBasicCellInfo();
	c.CalcCoulombEnergy();
	c.CalcCoulombDerivative();
	c.CalcLatticeDerivative();
	c.ShowEnergyDerivative();


	c.Finalise();
	return 0;
}
