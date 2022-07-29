#ifndef __CELL
#define __CELL

#include <Eigen/Dense>	// Matrix/Vector Arithematics?
#include <string>
#include <chrono>

#include "Atom.hpp"

#define DEF_MAX_ATOM_NUMBER 1024
#define DEF_PERIODIC_SUM_ACCURACY 10E-18
#define DEF_PERIODIC_SUM_WEIGHT   0.0123

class Cell
{

friend class Manager;	// Manager can access Cell privates

private:

const double TO_EV = 14.39964390675221758120;

Eigen::Vector3d lattice_param;
Eigen::Vector3d lattice_angle;

Eigen::Matrix3d lattice_matrix;		// lattice matrix {{a1,a2,a3},{b1,b2,b3},{c1,c2,c3}}

Eigen::Vector3d real_vector[3];		// Accessing example : real_vector[i](0,1,2)
Eigen::Vector3d reci_vector[3];

double volume;				// Cell Volume

int NumberOfAtoms;
Atom* AtomList[DEF_MAX_ATOM_NUMBER];	// New Del

// PERIODIC SUM PARAMETERS
double sigma;							// Sigma
const double accuracy = DEF_PERIODIC_SUM_ACCURACY;		// A parameter 
const double weight   = DEF_PERIODIC_SUM_WEIGHT;
double rcut,gcut;						// Cutoff tolerances
int h_max,k_max,l_max;
int ih_max,ik_max,il_max;


// Monopole Energy - i.e., point charge - e.g., core, shell, lone pair core ...
double mono_real_energy, mono_reci_energy, mono_reci_self_energy;
double mono_total_energy;

// Cell Strain Derivatives
Eigen::Matrix3d lattice_sd;	// ordering convention - e11(xx), e22(yy), e33(zz), e12(xy), e13(xz), e23(yz)
				//			 e1       e2       e3       e6       e5       e4
// PERIODIC SUM WORKING VARIABLES


// Wtime
std::chrono::duration<double> energy_real_wtime;
std::chrono::duration<double> energy_reci_wtime;
std::chrono::duration<double> derivative_real_wtime;		// std::cout << this->energy_wtime.count() << " s\n";
std::chrono::duration<double> derivative_reci_wtime;		// std::cout << this->energy_wtime.count() << " s\n";
// Lattice Sum Info
int energy_real_sum_cnt;	// Number of Operation Counter [1] - Energy, [2] - Derivatives
int energy_reci_sum_cnt;	// Number of Operation Counter [1] - Energy, [2] - Derivatives
int derivative_real_sum_cnt;
int derivative_reci_sum_cnt;

public:

Cell( std::string );

void CalcCoulombEnergy();					// () field can potentially be used for adding constraints, e.g., E fields later by overloading
void CalcCoulombDerivative();

void ShowBasicCellInfo() const;
void ShowEnergyDerivative() const;

virtual ~Cell();

};

#endif
