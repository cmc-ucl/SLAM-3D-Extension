#ifndef __LONEPAIR_MATRIX_H
#define __LONEPAIR_MATRIX_H

#define EV_UNIT (14.39964390675221758120)
#define TO_BOHR_RADII (0.52917721067)							// div Angstrom -> Bohr, mul Bohr -> Angstrom
#define HA_TO_EV_UNIT (EV_UNIT/TO_BOHR_RADII)						// mul Ha       -> eV
#define FHA_TO_FEV_UNIT (EV_UNIT/TO_BOHR_RADII/TO_BOHR_RADII)				// mul Ha/Bohr  -> eV/Angstrom
#define FFHA_TO_FFEV_UNIT (EV_UNIT/TO_BOHR_RADII/TO_BOHR_RADII/TO_BOHR_RADII)		// mul Ha/Bohr^2-> eV/Angstrom^2

#include <Eigen/Core>
#include "Integral_lib.hpp"

class LonePairMatrix
{
public:

	Eigen::Matrix4d transform_matrix;
	const int b_serach( const double dist, const std::vector<double>& integral_knot );
	const Eigen::Matrix4d& GetTransformationMatrix( const Eigen::Vector3d& Rij );
};

class LonePairMatrix_H : public LonePairMatrix	// 'public' specification is required ... take the public features as public
{

public: 

///	///	///	///	///	///	///	///
	void test()				// Binary Link Test
	{	printf("LonePairMatrix@@@@\n");
	}
	void test2();

	// Test Function for Validation
	const double NIntegral_test_real( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4] );
///	///	///	///	///	///	///	///
};

#endif
