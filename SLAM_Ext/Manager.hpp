#ifndef __INTERACTION_Manager
#define __INTERACTION_Manager

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

#include "Cell.hpp"
//#include "Atom.hpp"		// Exiplicit declaration
#include "LonePairMatrix.hpp"

class Manager // Interaction Manager
{

std::vector<double> man_vec;
LonePairMatrix_H man_lp_matrix_h;

public:

void InitialiseEnergy( Cell& C );
void InitialiseDerivative( Cell& C );						// Eigen::Vec/Mat::setZero(); (for addition assigment operations)
void InitialisePeriodicSysParameter( Cell& C );					// Parameters ... Spliting Workload of Lattice Sum (Real/Reciprocal)	// Under Dev 28 Jul Thurs 2022

// Monopole - Monopole Interaction (charge vs charge)
// decltype intrinsic - tells 'core-core', 'shell-core', 'core-shgl'
void CoulombMonoMonoReal( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector );
void CoulombMonoMonoSelf( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector );
void CoulombMonoMonoReci( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector );

// First Geometric Derivatives
void CoulombDerivativeReal( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector );
void CoulombDerivativeSelf( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector );
void CoulombDerivativeReci( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector );

// Strain Derivatives
void StrainDerivativeReal( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector );
void StrainDerivativeSelf( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector );
void StrainDerivativeReci( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector );

//
void InitialiseLonePairEnergy( Cell& C );
void InitialiseLonePairDerivative( Cell& C );

void InitialiseSCF();
bool IsSCFDone( const double tol );

void GetLonePairGroundState( Cell& C );

/* Return Transformation Matrix - Transformation Matrix (T) rotates global symmetry to local, e.g., Hij = Tai Tbj Hab'	- inverse transform of H' (local) matrix */
//const Eigen::Matrix4d& LonePairGetTransformationMatrix( Eigen::Matrix4d& transform_matrix /*in/out*/, const Eigen::Vector3d cart_i, const Eigen::Vector3d cart_j );

// Interaction - Lone Pair
void CoulombLonePairReal( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector, const bool is_first_scf );	// scf flag for calculating 
void CoulombLonePairSelf( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector, const bool is_first_scf );
void CoulombLonePairReci( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector, const bool is_first_scf );

void CoulombLonePairDerivativeReal( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector );
void CoulombLonePairDerivativeSelf( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector );
void CoulombLonePairDerivativeReci( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector );

void StrainLonePairDerivativeReal( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector );
void StrainLonePairDerivativeSelf( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector );
void StrainLonePairDerivativeReci( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector );

};

#endif
