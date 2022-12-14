#ifndef __CELL
#define __CELL

#include <Eigen/Dense>
#include <string>
#include <chrono>
#include <omp.h>

#include "Atom.hpp"

#define DEF_MAX_ATOM_NUMBER 1024
#define DEF_PERIODIC_SUM_ACCURACY 10E-18

//#define DEF_PERIODIC_SUM_WEIGHT   1.234
//#define DEF_PERIODIC_SUM_WEIGHT   5.678
//#define DEF_PERIODIC_SUM_WEIGHT   0.0123

// Benchmark done ... doesn't lose accuracy in real/reciprocal calculations
//#define DEF_PERIODIC_SUM_WEIGHT   1.123
//#define DEF_PERIODIC_SUM_WEIGHT   1.00141	// similar to gulp core only - rspeed 1.123^3

#define DEF_PERIODIC_SUM_WEIGHT   1.123

#define LONEPAIR_SCF_MAX 50
#define LONEPAIR_SCF_TOL 10E-9
#define LONEPAIR_DERIVATIVE_TOL 10E-8

class Cell
{

friend class Manager;	// Manager can access Cell privates

private:

const double TO_EV = 14.39964390675221758120;

Eigen::Vector3d lattice_param;
Eigen::Vector3d lattice_angle;

Eigen::Matrix3d lattice_matrix;		// lattice matrix {{a1,a2,a3},{b1,b2,b3},{c1,c2,c3}}
Eigen::Matrix3d lattice_matrix_transpose;
Eigen::Matrix3d lattice_matrix_inverse;

Eigen::Vector3d real_vector[3];		// Accessing example : real_vector[i](0,1,2)
Eigen::Vector3d reci_vector[3];

// Cell Strain Derivatives
Eigen::Matrix3d lattice_sd;	// ordering convention - e11(xx), e22(yy), e33(zz), e12(xy), e13(xz), e23(yz) - e1, e2, e3, e6, e5, e4
// Cell Lattice Derivatives
Eigen::Matrix3d lattice_derivative;	// ordering convention - Not Determined yet

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

/*
// LP PERIODIC SUM PARAMETERS
double lp_sigma;						// Sigma
double lp_rcut,lp_gcut;						// Cutoff tolerances
int lp_h_max,lp_k_max,lp_l_max;
int lp_ih_max,lp_ik_max,lp_il_max;
*/

// Monopole Energy - i.e., point charge - e.g., core, shell, lone pair core ...
double mono_real_energy, mono_reci_energy, mono_reci_self_energy;
double mono_total_energy;

// LonePair Energy Contributions ... Remember some terms below only contain the half of the contributions
double lp_real_energy, lp_reci_energy;
double lp_eval_sum;
double lp_total_energy; // real + reci + evals

// LONEPAIR__ Calculatiosn LonePiar involved
double lp_scf_energy_prev;				// SCF usage ... logging prev scf energy of lone pairs
double lp_scf_energy_curr;				// SCF usage ... logging current scf energy of lone pairs
const int scf_iter_max = LONEPAIR_SCF_MAX;


// Wtime Variables - MM calculation
std::chrono::duration<double> energy_real_wtime;
std::chrono::duration<double> energy_reci_wtime;
std::chrono::duration<double> derivative_real_wtime;		// std::cout << this->energy_wtime.count() << " s\n";
std::chrono::duration<double> derivative_reci_wtime;		// std::cout << this->energy_wtime.count() << " s\n";
// Lattice Sum Info - MM calculation
int energy_real_sum_cnt;	// Number of Operation Counter [1] - Energy, [2] - Derivatives
int energy_reci_sum_cnt;	// Number of Operation Counter [1] - Energy, [2] - Derivatives
int derivative_real_sum_cnt;
int derivative_reci_sum_cnt;

// Wtime Variables - LP calculation
std::chrono::duration<double> lp_energy_real_wtime;
std::chrono::duration<double> lp_energy_reci_wtime;
std::chrono::duration<double> lp_derivative_real_wtime;
std::chrono::duration<double> lp_derivative_reci_wtime;

public:

Cell( std::string );					// Initialiser

// Calculations - Shell/Core classics?
virtual void CalcCoulombEnergy();			// () field can potentially be used for adding constraints, e.g., E fields later by overloading
virtual void CalcCoulombDerivative();
virtual void CalcLatticeDerivative();

// LONEPAIR__ Calculatiosn LonePiar involved
virtual void CalcLonePairCoulombEnergy();
virtual void CalcLonePairCoulombDerivative();

// Print things
virtual void ShowBasicCellInfo() const;
virtual void ShowEnergyDerivative() const;
virtual void Finalise() const;

// Destructor
virtual ~Cell();

};

#endif
