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

#define MX_C 32

class Manager // Interaction Manager
{

std::vector<double> man_vec;
LonePairMatrix_H man_lp_matrix_h;

// lp - pc (point charge)
Eigen::Matrix4d real_lp_h_pc[MX_C][MX_C];

Eigen::Matrix4d real_lp_h_pc_x[MX_C][MX_C];
Eigen::Matrix4d real_lp_h_pc_y[MX_C][MX_C];
Eigen::Matrix4d real_lp_h_pc_z[MX_C][MX_C];

// lp - ps (point shell)
Eigen::Matrix4d real_lp_h_ps[2][MX_C][MX_C];

Eigen::Matrix4d real_lp_h_ps_x[2][MX_C][MX_C];
Eigen::Matrix4d real_lp_h_ps_y[2][MX_C][MX_C];
Eigen::Matrix4d real_lp_h_ps_z[2][MX_C][MX_C];

// lp - lp (core)
Eigen::Matrix4d real_lp_h_lc[MX_C][MX_C];

Eigen::Matrix4d real_lp_h_lc_x[MX_C][MX_C];
Eigen::Matrix4d real_lp_h_lc_y[MX_C][MX_C];
Eigen::Matrix4d real_lp_h_lc_z[MX_C][MX_C];

// lp - lp
Eigen::Matrix4d real_lp_h_lp[MX_C];

Eigen::Matrix4d real_lp_h_lp_x[MX_C][MX_C];
Eigen::Matrix4d real_lp_h_lp_y[MX_C][MX_C];
Eigen::Matrix4d real_lp_h_lp_z[MX_C][MX_C];

Eigen::Matrix4d real_lp_h_lp_xx[MX_C][MX_C];
Eigen::Matrix4d real_lp_h_lp_xy[MX_C][MX_C];
Eigen::Matrix4d real_lp_h_lp_xz[MX_C][MX_C];
Eigen::Matrix4d real_lp_h_lp_yx[MX_C][MX_C];
Eigen::Matrix4d real_lp_h_lp_yy[MX_C][MX_C];
Eigen::Matrix4d real_lp_h_lp_yz[MX_C][MX_C];
Eigen::Matrix4d real_lp_h_lp_zx[MX_C][MX_C];
Eigen::Matrix4d real_lp_h_lp_zy[MX_C][MX_C];
Eigen::Matrix4d real_lp_h_lp_zz[MX_C][MX_C];

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

// Setting LP H matrices
const Eigen::Matrix4d& set_h_matrix_real_pc( LonePair* lp, const Eigen::Vector3d& R, const double sig, const int lp_i, const int pc_i );
void set_h_matrix_real_pc_derivative( LonePair* lp, const Eigen::Vector3d& R, const double sig, const int lp_i, const int pc_i );

//

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
