#include <iostream>
#include "cell.hpp"

#define TARGET_INPUT "geo.txt"

int main()
{
	Cell c(TARGET_INPUT);
	c.ShowBasicCellInfo();
	c.CalcCoulombEnergy();
	c.CalcCoulombDerivative();
	c.ShowEnergyDerivative();
	return 0;
}
