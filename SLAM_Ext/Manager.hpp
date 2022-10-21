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
LonePairMatrix_H man_lp_matrix_h;	// LonePairMatrix_H : Instance ... managing solving integration tools


////	////	////	////	////	////	////	////
// lp - pc (point charge)
Eigen::Matrix4d real_lp_h_pc[MX_C][MX_C];

Eigen::Matrix4d reci_lp_h_pc[MX_C][MX_C];

Eigen::Matrix4d real_lp_h_pc_x[MX_C][MX_C];
Eigen::Matrix4d real_lp_h_pc_y[MX_C][MX_C];
Eigen::Matrix4d real_lp_h_pc_z[MX_C][MX_C];

Eigen::Matrix4d reci_lp_h_pc_x[MX_C][MX_C];
Eigen::Matrix4d reci_lp_h_pc_y[MX_C][MX_C];
Eigen::Matrix4d reci_lp_h_pc_z[MX_C][MX_C];

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
////	////	////	////	////	////	////	////



Eigen::Matrix4d man_matrix4d_h_real_ws[2];			// Workspace H matrix real part
Eigen::Matrix4d man_matrix4d_h_real_derivative_ws[3];		// Workspace H matrix real derivative
Eigen::Matrix4d man_matrix4d_h_real_derivative_out[3];
Eigen::Matrix4d man_matrix4d_h_real_derivative2_ws[9];		// Workspace H matrix real derivative2
Eigen::Matrix4d man_matrix4d_h_real_derivative2_out[9];


/* 
	Storage for H Matrix (E) - LP vs PointCharge Interaction (core/shell/lp core)
	Saving Interaction of a lone pair density (in the central sublattice) w.r.t. surrounding classical entities 
*/
Eigen::Matrix4d LPC_H_Real[MX_C][MX_C][2];	// [i][j][0] ... if 'j' is core or lp_core // [i][j][1] ... if 'j' is shell
Eigen::Matrix4d LPC_H_Reci[MX_C][MX_C][2];	// same convention ... if 'i=j' with h'=k'=l' - reciprocal self (i.e., lone pair electron interacting with its core in the central sublattice)

// Storage for H Matrix (E) - LP vs LP
Eigen::Matrix4d LPLP_H_Real[MX_C][MX_C];	// Save Monopolar Contribution LP <---> LP
Eigen::Matrix4d LPLP_H_Real_x[MX_C][MX_C];	// Save Dipolar   Contribution LP <---> LP
Eigen::Matrix4d LPLP_H_Real_y[MX_C][MX_C];
Eigen::Matrix4d LPLP_H_Real_z[MX_C][MX_C];
Eigen::Matrix4d LPLP_H_Reci[MX_C][MX_C];


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

// RealSpace - Setting LP H Matrices
const Eigen::Matrix4d& set_h_matrix_real_pc( LonePair* lp, const Eigen::Vector3d& R, const double sig, const int lp_i, const int pc_i );
void set_h_matrix_real_pc_derivative( LonePair* lp, const Eigen::Vector3d& R, const double sig, const int lp_i, const int pc_i );
void set_h_matrix_real_pc_derivative2( LonePair* lp, const Eigen::Vector3d& R, const double sig, const int lp_i, const int pc_i );

// ReciSpace - Setting LP H Matrices CosinePart
const Eigen::Matrix4d& set_h_matrix_reci_cos( LonePair* lp, const Eigen::Vector3d& G, const int lp_i, const int pc_i );
void set_h_matrix_reci_derivative_cos( LonePair* lp, const Eigen::Vector3d& G, const int lp_i, const int pc_i );

// ReciSpace - Setting LP H Matrices SinePart
const Eigen::Matrix4d& set_h_matrix_reci_sin( LonePair* lp, const Eigen::Vector3d& G, const int lp_i, const int pc_i );
void set_h_matrix_reci_derivative_sin( LonePair* lp, const Eigen::Vector3d& G, const int lp_i, const int pc_i );



// Internal Uses in " CoulombLonePairReal " below
void support_h_matrix_real( const LonePair* lp, const double& sigma, const Eigen::Vector3d& Rij, /* in/out */ Eigen::Matrix4d& h_mat_ws, Eigen::Matrix4d& h_mat_out );
void support_h_matrix_real_derivative( const LonePair* lp, const double& sigma, const Eigen::Vector3d& Rij, /* workspace */ Eigen::Matrix4d (&h_mat_ws)[3], /* out */ Eigen::Matrix4d (&h_mat_out)[3] );
// ACTUAL USE - H (E) Matrix Calculations
void set_h_matrix_real( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector );







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
