#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include <omp.h>
#define NUM_OMP_THREAD 4

#include "LonePairMatrix.hpp"

// Internal Use Only
double radial( double a, double b, double c, double d, double r ) { return a*r*r*r + b*r*r + c*r + d; }
double grid( double delta_r ) { return ceil(fabs(GRID_W/log(delta_r))); }
//

int LonePairMatrix::b_search( const double dist, const std::vector<double>& integral_knot )
{
	const int knot_stride = integral_knot.size();

    int le = 0;
    int re = knot_stride - 1;

    while(1)
    {   if( integral_knot[(le+re)/2] < dist )
            le = (le+re)/2;
        else // dist <= integral_knot[(le+re)/2]
            re = (le+re)/2;

        if( (le+1) == re )
            break;
    }

    return le;		// Return knot index 
}

const Eigen::Matrix4d& LonePairMatrix::GetTransformationMatrix( const Eigen::Vector3d& Rij )
//const Eigen::Matrix4d& LonePairMatrix::GetTransformationMatrix( Eigen::Vector3d& Rij )
{                                                               // reserved for : Cell.lp_transformation_matrix, AtomList[i/j]->cart, static_cast<Shell*>(AtomList[i/j])->shel_cart,
        //Eigen::Vector3d Rij = cart_i - cart_j;
        //Eigen::Vector3d Rij = cart_j - cart_i;                  // Using This Convention ... Must follow that 'cart_i' has to be a core position of LonePair of interest
                                                                //                                            'cart_j' has to be a classic core / or the other LP
	//Eigen::Vector3d Rij = -r;	
	//Eigen::Vector3d Rij = r;	

        const double Rxy = sqrt(Rij(0)*Rij(0) + Rij(1)*Rij(1));
        const double R   = Rij.norm();
        double n1, n2, tmp;// dummy variables for workspace

        this->transform_matrix.setZero();// init Return

        this->transform_matrix(0,0) = 1.;

        //if( (Rij(0) == 0) && (Rij(1) == 0) && (Rij(2) > 0) )                    // if vector 'Rij' is on z-axis
        if( ( fabs(Rij(0)) < 10E-11 ) && ( fabs(Rij(1)) < 10E-11 ) && (Rij(2) > 0) )                    // if vector 'Rij' is on z-axis
        {       this->transform_matrix(1,1) = 1.;
                this->transform_matrix(2,2) = 1.;
                this->transform_matrix(3,3) = 1.;                             // set lower-right block 3x3 matrix as I
        }
        //else if( (Rij(0) == 0) && (Rij(1) == 0) && (Rij(2) < 0) )               // if vector 'Rij' is on negative z-axis
        else if( ( fabs(Rij(0)) < 10E-11 ) && ( fabs(Rij(1)) < 10E-11 ) && (Rij(2) < 0) )               // if vector 'Rij' is on negative z-axis
        {       this->transform_matrix(1,1) = 1.;
                this->transform_matrix(2,2) = 1.;
                this->transform_matrix(3,3) =-1.;                             // set the matrix has xy-plane reflection
        }
        else
        {       this->transform_matrix(3,1) = Rij(0)/R;
                this->transform_matrix(3,2) = Rij(1)/R;
                this->transform_matrix(3,3) = Rij(2)/R;                       // set local z-axis (k') in the transformed symmetry

                this->transform_matrix(2,1) = Rij(2)*Rij(0)/Rxy;
                this->transform_matrix(2,2) = Rij(2)*Rij(1)/Rxy;
                this->transform_matrix(2,3) = -R*sqrt(1.-Rij(2)*Rij(2)/R/R);

                n1 = 1./sqrt(this->transform_matrix(2,1)*this->transform_matrix(2,1) + this->transform_matrix(2,2)*this->transform_matrix(2,2) + this->transform_matrix(2,3)*this->transform_matrix(2,3));

                for(int k=0;k<3;k++)
                {       tmp = this->transform_matrix(2,k+1);
                        tmp = tmp*n1;
                        this->transform_matrix(2,k+1) = tmp;                  // set local y-axis (j') in the transformed symmetry
                }

                this->transform_matrix(1,1) = n1/R * ( Rij(2)*Rij(2)*Rij(1)/Rxy + R*Rij(1)*sqrt(1.-Rij(2)*Rij(2)/R/R) );
                this->transform_matrix(1,2) = n1/R * (-Rij(2)*Rij(2)*Rij(0)/Rxy - R*Rij(0)*sqrt(1.-Rij(2)*Rij(2)/R/R) );
                this->transform_matrix(1,3) = 0.;

                n2 = 1./sqrt(this->transform_matrix(1,1)*this->transform_matrix(1,1) + this->transform_matrix(1,2)*this->transform_matrix(1,2) + this->transform_matrix(1,3)*this->transform_matrix(1,3));

                for(int k=0;k<3;k++)
                {       tmp = this->transform_matrix(1,k+1);
                        tmp = tmp/n2;
                        this->transform_matrix(1,k+1) = tmp;
                }
        }

	// Saving the LowerBlock Diagonal 3x3 Matrix
	for(int i=0;i<3;i++)
	{	for(int j=0;j<3;j++) { this->transform_matrix_shorthand(i,j) = this->transform_matrix(i+1,j+1); }}
	// Return
        return this->transform_matrix;
}

////	////	////	////	////	////

////	Real Space Position Integral

///	////	////	////	////	////

double LonePairMatrix_H::real_position_integral( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4] )
{
	double res = 0.;

	for(int i=0;i<integral_knot.size()-1;i++)
	{
		res += (-0.5*(Rp[3][i]*Rs[3][i]*pow(integral_knot[i],2)) - ((Rs[2][i]*Rp[3][i] + Rp[2][i]*Rs[3][i])*pow(integral_knot[i],3))/3. - ((Rp[2][i]*Rs[2][i]
		+ Rs[1][i]*Rp[3][i] + Rp[1][i]*Rs[3][i])*pow(integral_knot[i],4))/4. - ((Rs[1][i]*Rp[2][i] + Rp[1][i]*Rs[2][i] + Rs[0][i]*Rp[3][i]
		+ Rp[0][i]*Rs[3][i])*pow(integral_knot[i],5))/5. - ((Rp[1][i]*Rs[1][i] + Rs[0][i]*Rp[2][i] + Rp[0][i]*Rs[2][i])*pow(integral_knot[i],6))/6.
		- ((Rs[0][i]*Rp[1][i] + Rp[0][i]*Rs[1][i])*pow(integral_knot[i],7))/7. - (Rp[0][i]*Rs[0][i]*pow(integral_knot[i],8))/8. + (Rp[3][i]*Rs[3][i]*pow(integral_knot[i+1],2))/2.
		+ ((Rs[2][i]*Rp[3][i] + Rp[2][i]*Rs[3][i])*pow(integral_knot[i+1],3))/3. + ((Rp[2][i]*Rs[2][i] + Rs[1][i]*Rp[3][i] + Rp[1][i]*Rs[3][i])*pow(integral_knot[i+1],4))/4.
		+ ((Rs[1][i]*Rp[2][i] + Rp[1][i]*Rs[2][i] + Rs[0][i]*Rp[3][i] + Rp[0][i]*Rs[3][i])*pow(integral_knot[i+1],5))/5. + ((Rp[1][i]*Rs[1][i] + Rs[0][i]*Rp[2][i]
		+ Rp[0][i]*Rs[2][i])*pow(integral_knot[i+1],6))/6. + ((Rs[0][i]*Rp[1][i] + Rp[0][i]*Rs[1][i])*pow(integral_knot[i+1],7))/7. + (Rp[0][i]*Rs[0][i]*pow(integral_knot[i+1],8))/8.)/sqrt(3);
	}
	
	// distance unit
	return res*TO_BOHR_RADII;
}

////	////	////	////	////	////

////	Reci Space Self Integral : Applies for LPe - LPcore self interaction when 'i'LPe / 'j' LPcore i = j; case looping through the central sublattice

///	////	////	////	////	////

double LonePairMatrix_H::reci_self_integral_ss( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;	

	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);

		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
		for(int k=0;k<(int)mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*erf(r/sig)/r;
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*erf(r/sig)/r;
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res*HA_TO_EV_UNIT;
}

// xx = yy = zz
double LonePairMatrix_H::reci_self_integral_xx( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;	
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);

		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
		for(int k=0;k<(int)mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*erf(r/sig)/r;
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*erf(r/sig)/r;
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res*HA_TO_EV_UNIT;
}

////	////	////	////	////	////

////	Reci Space Self Integral d/dRi ... (derivative w.r.t. 'i' core of a LP)

////	grad x - sx = grad y - sy = grad z - sz

////	////	////	////	////	////
//((2*r)/(pow(M_E,pow(r,2)/pow(sig,2))*sqrt(M_PI)*sig) - erf(r/sig))/(sqrt(3)*pow(r,2))
double LonePairMatrix_H::reci_self_integral_sx_grad_x( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);

		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
		for(int k=0;k<(int)mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*(((2*r)/(pow(M_E,pow(r,2)/pow(sig,2))*sqrt(M_PI)*sig) - erf(r/sig))/(sqrt(3)*pow(r,2)));
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*(((2*r)/(pow(M_E,pow(r,2)/pow(sig,2))*sqrt(M_PI)*sig) - erf(r/sig))/(sqrt(3)*pow(r,2)));
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res*HA_TO_EV_UNIT; // Ha / a
}



////	////	////	////	////	////

////	Real Space Integrals

////	////	////	////	////	////

double LonePairMatrix_H::real_ss_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;

	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);

		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
		for(int k=0;k<(int)mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_ss(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_ss(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res*HA_TO_EV_UNIT;
}

double LonePairMatrix_H::real_sz_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
		for(int k=0;k<(int)mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_sz(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_sz(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res*HA_TO_EV_UNIT;
}

double LonePairMatrix_H::real_xx_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
		for(int k=0;k<(int)mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_xx(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_xx(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res*HA_TO_EV_UNIT;
}

double LonePairMatrix_H::real_zz_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;

	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
		for(int k=0;k<(int)mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_zz(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_zz(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res*HA_TO_EV_UNIT;
}

double LonePairMatrix_H::real_sx_grad_x_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
		for(int k=0;k<(int)mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_x_sx(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_x_sx(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res*FHA_TO_FEV_UNIT;
}

double LonePairMatrix_H::real_xz_grad_x_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4],  double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
		for(int k=0;k<(int)mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_x_xz(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_x_xz(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res*FHA_TO_FEV_UNIT;
}

//double LonePairMatrix_H::real_ss_grad_z_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], const double sig, double d )
double LonePairMatrix_H::real_ss_grad_z_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;

	const int knot_d = LonePairMatrix::b_search( d, integral_knot );	// Get 'd' location on a knot

	for(int i=0;i<integral_knot.size()-1;i++)
	{
		if( i != knot_d )
		{
			dr   = integral_knot[i+1] - integral_knot[i];
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);
	
			#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
			for(int k=0;k<(int)mesh;k++)
			{
				r  = integral_knot[i] + k * r_inc;
				fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_derivative_z_ss(sig,r,d);
				r  = integral_knot[i] + (k+1) * r_inc;
				fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_derivative_z_ss(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
		}
		else
		{	
			dr   = d - integral_knot[i];
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);

			#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
			for(int k=0;k<(int)mesh;k++)
			{
				r  = integral_knot[i] + k * r_inc;
				fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_derivative_z_ss_left(sig,r,d);
				r  = integral_knot[i] + (k+1) * r_inc;
				fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_derivative_z_ss_left(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
			
			dr   = integral_knot[i+1] - d;
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);

			#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
			for(int k=0;k<(int)mesh;k++)
			{
				r  = d + k * r_inc;
				fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_derivative_z_ss_right(sig,r,d);
				r  = d + (k+1) * r_inc;
				fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_derivative_z_ss_right(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
		}
	}
	// eV Unit
	if( std::isnan(res) ){	printf("ss grad z - isnan\n"); exit(1); }

	return res*FHA_TO_FEV_UNIT;
}

double LonePairMatrix_H::real_sz_grad_z_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;

	const int knot_d = LonePairMatrix::b_search( d, integral_knot );	// Get 'd' location on a knot

	for(int i=0;i<integral_knot.size()-1;i++)
	{
		if( i != knot_d )
		{
			dr   = integral_knot[i+1] - integral_knot[i];
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);
		
			#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
			for(int k=0;k<(int)mesh;k++)
			{
				r  = integral_knot[i] + k * r_inc;
				fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_z_sz(sig,r,d);
				r  = integral_knot[i] + (k+1) * r_inc;
				fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_z_sz(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
		}
		else
		{	
			dr   = d - integral_knot[i];
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);

			#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
			for(int k=0;k<(int)mesh;k++)
			{
				r  = integral_knot[i] + k * r_inc;
				fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_z_sz_left(sig,r,d);
				r  = integral_knot[i] + (k+1) * r_inc;
				fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_z_sz_left(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
			
			dr   = integral_knot[i+1] - d;
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);

			#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
			for(int k=0;k<(int)mesh;k++)
			{
				r  = d + k * r_inc;
				fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_z_sz_right(sig,r,d);
				r  = d + (k+1) * r_inc;
				fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_z_sz_right(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
		}
	}
	if( std::isnan(res) ){	printf("sz grad z - isnan\n"); exit(1); }
	// eV Unit
	return res*FHA_TO_FEV_UNIT;
}

double LonePairMatrix_H::real_xx_grad_z_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
		for(int k=0;k<(int)mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_z_xx(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_z_xx(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res*FHA_TO_FEV_UNIT;
}


double LonePairMatrix_H::real_zz_grad_z_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;

	const int knot_d = LonePairMatrix::b_search( d, integral_knot );	// Get 'd' location on a knot

	for(int i=0;i<integral_knot.size()-1;i++)
	{
		if( i != knot_d )
		{
			dr   = integral_knot[i+1] - integral_knot[i];
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);
		
			#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
			for(int k=0;k<(int)mesh;k++)
			{
				r  = integral_knot[i] + k * r_inc;
				fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_z_zz(sig,r,d);
				r  = integral_knot[i] + (k+1) * r_inc;
				fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_z_zz(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
		}
		else
		{	
			dr   = d - integral_knot[i];
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);

			#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
			for(int k=0;k<(int)mesh;k++)
			{
				r  = integral_knot[i] + k * r_inc;
				fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_z_zz_left(sig,r,d);
				r  = integral_knot[i] + (k+1) * r_inc;
				fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_z_zz_left(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
			
			dr   = integral_knot[i+1] - d;
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);

			#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
			for(int k=0;k<(int)mesh;k++)
			{
				r  = d + k * r_inc;
				fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_z_zz_right(sig,r,d);
				r  = d + (k+1) * r_inc;
				fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative_z_zz_right(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
		}
	}
	if( std::isnan(res) ){	printf("zz grad z - isnan\n"); exit(1); }
	// eV Unit
	return res*FHA_TO_FEV_UNIT;
}

////	////	////	////	////	////

////	2nd derivatives w.r.t. point charges - RealSpace

////	////	////	////	////	////

//double et = -omp_get_wtime();
//et += omp_get_wtime();
//printf("calculated val : %20.12lf\n",res);
//printf("elapsed time   : %20.12lf\n",et);
//exit(1);

double LonePairMatrix_H::real_ss_grad2_xx_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		for(int k=0;k<mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_derivative2_xx_ss(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_derivative2_xx_ss(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}

	// eV Unit
	return res*FFHA_TO_FFEV_UNIT;
}

double LonePairMatrix_H::real_sz_grad2_xx_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		for(int k=0;k<mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_xx_sz(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_xx_sz(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}

	// eV Unit
	return res*FFHA_TO_FFEV_UNIT;
}

double LonePairMatrix_H::real_xx_grad2_xx_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		for(int k=0;k<mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_xx_xx(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_xx_xx(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res*FFHA_TO_FFEV_UNIT;
}

double LonePairMatrix_H::real_yy_grad2_xx_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		for(int k=0;k<mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_xx_yy(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_xx_yy(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res*FFHA_TO_FFEV_UNIT;
}

double LonePairMatrix_H::real_zz_grad2_xx_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		for(int k=0;k<mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_xx_zz(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_xx_zz(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res*FFHA_TO_FFEV_UNIT;
}

double LonePairMatrix_H::real_xy_grad2_xy_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		for(int k=0;k<mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_xy_xy(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_xy_xy(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res*FFHA_TO_FFEV_UNIT;
}

////	Derivative XZ

double LonePairMatrix_H::real_sx_grad2_xz_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		for(int k=0;k<mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_xz_sx(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_xz_sx(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res*FFHA_TO_FFEV_UNIT;
}

double LonePairMatrix_H::real_xz_grad2_xz_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		for(int k=0;k<mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_xz_xz(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_xz_xz(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res*FFHA_TO_FFEV_UNIT;
}

////	Derivative ZZ

double LonePairMatrix_H::real_ss_grad2_zz_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;

	const int knot_d = LonePairMatrix::b_search( d, integral_knot );	// Get 'd' location on a knot

	for(int i=0;i<integral_knot.size()-1;i++)
	{
		if( i != knot_d )
		{
			dr   = integral_knot[i+1] - integral_knot[i];
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);
		
			for(int k=0;k<mesh;k++)
			{
				r  = integral_knot[i] + k * r_inc;
				fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_ss(sig,r,d);
				r  = integral_knot[i] + (k+1) * r_inc;
				fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_ss(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
		}
		else
		{	
			dr   = d - integral_knot[i];
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);

			for(int k=0;k<mesh;k++)
			{
				r  = integral_knot[i] + k * r_inc;
				fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_ss_left(sig,r,d);
				r  = integral_knot[i] + (k+1) * r_inc;
				fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_ss_left(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
			
			dr   = integral_knot[i+1] - d;
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);

			for(int k=0;k<mesh;k++)
			{
				r  = d + k * r_inc;
				fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_ss_right(sig,r,d);
				r  = d + (k+1) * r_inc;
				fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_ss_right(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
			// Discontinuous Point AUx
			res += radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],d)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],d)*real_derivative2_aux_grad_z_ss(sig,d,d);
		}
	}

	// eV Unit
	return res*FFHA_TO_FFEV_UNIT;
}


double LonePairMatrix_H::real_sz_grad2_zz_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;

	const int knot_d = LonePairMatrix::b_search( d, integral_knot );	// Get 'd' location on a knot

	for(int i=0;i<integral_knot.size()-1;i++)
	{
		if( i != knot_d )
		{
			dr   = integral_knot[i+1] - integral_knot[i];
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);
		
			for(int k=0;k<mesh;k++)
			{
				r  = integral_knot[i] + k * r_inc;
				fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_sz(sig,r,d);
				r  = integral_knot[i] + (k+1) * r_inc;
				fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_sz(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
		}
		else
		{	
			dr   = d - integral_knot[i];
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);

			for(int k=0;k<mesh;k++)
			{
				r  = integral_knot[i] + k * r_inc;
				fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_sz_left(sig,r,d);
				r  = integral_knot[i] + (k+1) * r_inc;
				fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_sz_left(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
			
			dr   = integral_knot[i+1] - d;
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);

			for(int k=0;k<mesh;k++)
			{
				r  = d + k * r_inc;
				fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_sz_right(sig,r,d);
				r  = d + (k+1) * r_inc;
				fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_sz_right(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
			// Discontinuous Point AUx
			res += radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],d)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],d)*real_derivative2_aux_grad_z_sz(sig,d,d);
		}
	}

	// eV Unit
	return res*FFHA_TO_FFEV_UNIT;
}


double LonePairMatrix_H::real_xx_grad2_zz_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		for(int k=0;k<mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_xx(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_xx(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}

		// 
	}
	// eV Unit
	return res*FFHA_TO_FFEV_UNIT;
}

double LonePairMatrix_H::real_zz_grad2_zz_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double sig, double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	d = d/TO_BOHR_RADII;
	// Distance to Bohr
	sig = sig/TO_BOHR_RADII;

	const int knot_d = LonePairMatrix::b_search( d, integral_knot );	// Get 'd' location on a knot

	for(int i=0;i<integral_knot.size()-1;i++)
	{
		if( i != knot_d )
		{
			dr   = integral_knot[i+1] - integral_knot[i];
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);
		
			for(int k=0;k<mesh;k++)
			{
				r  = integral_knot[i] + k * r_inc;
				fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_zz(sig,r,d);
				r  = integral_knot[i] + (k+1) * r_inc;
				fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_zz(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
		}
		else
		{	
			dr   = d - integral_knot[i];
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);

			for(int k=0;k<mesh;k++)
			{
				r  = integral_knot[i] + k * r_inc;
				fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_zz_left(sig,r,d);
				r  = integral_knot[i] + (k+1) * r_inc;
				fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_zz_left(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
			
			dr   = integral_knot[i+1] - d;
			mesh = grid(dr);
			r_inc= dr/static_cast<double>(mesh);

			for(int k=0;k<mesh;k++)
			{
				r  = d + k * r_inc;
				fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_zz_right(sig,r,d);
				r  = d + (k+1) * r_inc;
				fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_derivative2_zz_zz_right(sig,r,d);
				res += r_inc*(fa+fb)/2.;
			}
			// Discontinuous Point AUx
			res += radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],d)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],d)*real_derivative2_aux_grad_z_zz(sig,d,d);
		}
	}

	// eV Unit
	return res*FFHA_TO_FFEV_UNIT;
}


////	////	////	////	////	////	////	////

////	Reciprocal Space Integral - LP

////	////	////	////	////	////	////	////

////	Cosine Part

//	Integral ss
double LonePairMatrix_H::reci_ss_cos( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double g )
{
        double res = 0.;
        double fa, fb, dr, mesh;
        double r,r_inc;
	// Unit 'g' Angs^-1 ... to Bohr^-1
	g = g*TO_BOHR_RADII;
        // Distance to Bohr

        for(int i=0;i<integral_knot.size()-1;i++)
        {
                dr   = integral_knot[i+1] - integral_knot[i];
                mesh = grid(dr);
                r_inc= dr/static_cast<double>(mesh);

		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
                for(int k=0;k<(int)mesh;k++)
                {
                        r  = integral_knot[i] + k * r_inc;
                        fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_reci_ss_cos(r,g);
                        r  = integral_knot[i] + (k+1) * r_inc;
                        fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_reci_ss_cos(r,g);
                        res += r_inc*(fa+fb)/2.;
                }

                // 
        }
        // Dimensionless ..
        return res;	
}

//	Integral xx <=> yy
double LonePairMatrix_H::reci_xx_cos( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double g )
{
        double res = 0.;
        double fa, fb, dr, mesh;
        double r,r_inc;
	// Unit 'g' Angs^-1 ... to Bohr^-1
	g = g*TO_BOHR_RADII;
        // Distance to Bohr

        for(int i=0;i<integral_knot.size()-1;i++)
        {
                dr   = integral_knot[i+1] - integral_knot[i];
                mesh = grid(dr);
                r_inc= dr/static_cast<double>(mesh);

		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
                for(int k=0;k<(int)mesh;k++)
                {
                        r  = integral_knot[i] + k * r_inc;
                        fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_reci_xx_cos(r,g);
                        r  = integral_knot[i] + (k+1) * r_inc;
                        fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_reci_xx_cos(r,g);
                        res += r_inc*(fa+fb)/2.;
                }

                // 
        }
        // Dimensionless
        return res;
}

//	Integral zz
double LonePairMatrix_H::reci_zz_cos( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double g )
{
        double res = 0.;
        double fa, fb, dr, mesh;
        double r,r_inc;
	// Unit 'g' Angs^-1 ... to Bohr^-1
	g = g*TO_BOHR_RADII;
        // Distance to Bohr

        for(int i=0;i<integral_knot.size()-1;i++)
        {
                dr   = integral_knot[i+1] - integral_knot[i];
                mesh = grid(dr);
                r_inc= dr/static_cast<double>(mesh);

		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
                for(int k=0;k<(int)mesh;k++)
                {
                        r  = integral_knot[i] + k * r_inc;
                        fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_reci_zz_cos(r,g);
                        r  = integral_knot[i] + (k+1) * r_inc;
                        fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_reci_zz_cos(r,g);
                        res += r_inc*(fa+fb)/2.;
                }

                // 
        }
        // Dimensionless
        return res;
}

////	Sine Part

//	Integral sz
double LonePairMatrix_H::reci_sz_sin( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double g )
{
        double res = 0.;
        double fa, fb, dr, mesh;
        double r,r_inc;
	// Unit 'g' Angs^-1 ... to Bohr^-1
	g = g*TO_BOHR_RADII;
        // Distance to Bohr

        for(int i=0;i<integral_knot.size()-1;i++)
        {
                dr   = integral_knot[i+1] - integral_knot[i];
                mesh = grid(dr);
                r_inc= dr/static_cast<double>(mesh);

		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
                for(int k=0;k<(int)mesh;k++)
                {
                        r  = integral_knot[i] + k * r_inc;
                        fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_reci_sz_sin(r,g);
                        r  = integral_knot[i] + (k+1) * r_inc;
                        fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_reci_sz_sin(r,g);
                        res += r_inc*(fa+fb)/2.;
                }

                // 
        }
        // Dimensionledd
        return res;
}

////	////	////

//// 	Reci Space Integral - LP derivative w.r.t. 'g' vector components

////	////	////

////	Cosine Part

//	Integral xz 'd/dgx'
double LonePairMatrix_H::reci_xz_grad_gx_cos( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double g )
{
        double res = 0.;
        double fa, fb, dr, mesh;
        double r,r_inc;
	// Unit 'g' Angs^-1 ... to Bohr^-1
	g = g*TO_BOHR_RADII;
        // Distance to Bohr

        for(int i=0;i<integral_knot.size()-1;i++)
        {
                dr   = integral_knot[i+1] - integral_knot[i];
                mesh = grid(dr);
                r_inc= dr/static_cast<double>(mesh);

		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
                for(int k=0;k<(int)mesh;k++)
                {
                        r  = integral_knot[i] + k * r_inc;
                        fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_reci_derivative_x_xz_cos(r,g);
                        r  = integral_knot[i] + (k+1) * r_inc;
                        fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_reci_derivative_x_xz_cos(r,g);
                        res += r_inc*(fa+fb)/2.;
                }

                // 
        }
        // Bohr -> Angstrom
        return res*TO_BOHR_RADII;
}

//	Integral ss 'd/dgz'
double LonePairMatrix_H::reci_ss_grad_gz_cos( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double g )
{
        double res = 0.;
        double fa, fb, dr, mesh;
        double r,r_inc;
	// Unit 'g' Angs^-1 ... to Bohr^-1
	g = g*TO_BOHR_RADII;
        // Distance to Bohr

        for(int i=0;i<integral_knot.size()-1;i++)
        {
                dr   = integral_knot[i+1] - integral_knot[i];
                mesh = grid(dr);
                r_inc= dr/static_cast<double>(mesh);

		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
                for(int k=0;k<(int)mesh;k++)
                {
                        r  = integral_knot[i] + k * r_inc;
                        fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_reci_derivative_z_ss_cos(r,g);
                        r  = integral_knot[i] + (k+1) * r_inc;
                        fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_reci_derivative_z_ss_cos(r,g);
                        res += r_inc*(fa+fb)/2.;
                }

                // 
        }
        // Bohr -> Angstrom
        return res*TO_BOHR_RADII;
}

//	Integral xx 'd/dgz'
double LonePairMatrix_H::reci_xx_grad_gz_cos( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double g )
{
        double res = 0.;
        double fa, fb, dr, mesh;
        double r,r_inc;
	// Unit 'g' Angs^-1 ... to Bohr^-1
	g = g*TO_BOHR_RADII;
        // Distance to Bohr

        for(int i=0;i<integral_knot.size()-1;i++)
        {
                dr   = integral_knot[i+1] - integral_knot[i];
                mesh = grid(dr);
                r_inc= dr/static_cast<double>(mesh);

		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
                for(int k=0;k<(int)mesh;k++)
                {
                        r  = integral_knot[i] + k * r_inc;
                        fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_reci_derivative_z_xx_cos(r,g);
                        r  = integral_knot[i] + (k+1) * r_inc;
                        fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_reci_derivative_z_xx_cos(r,g);
                        res += r_inc*(fa+fb)/2.;
                }

                // 
        }
        // Bohr -> Angstrom
        return res*TO_BOHR_RADII;
}

//	Integral zz 'd/dgz'
double LonePairMatrix_H::reci_zz_grad_gz_cos( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double g )
{
        double res = 0.;
        double fa, fb, dr, mesh;
        double r,r_inc;
	// Unit 'g' Angs^-1 ... to Bohr^-1
	g = g*TO_BOHR_RADII;
        // Distance to Bohr

        for(int i=0;i<integral_knot.size()-1;i++)
        {
                dr   = integral_knot[i+1] - integral_knot[i];
                mesh = grid(dr);
                r_inc= dr/static_cast<double>(mesh);

		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
                for(int k=0;k<(int)mesh;k++)
                {
                        r  = integral_knot[i] + k * r_inc;
                        fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_reci_derivative_z_zz_cos(r,g);
                        r  = integral_knot[i] + (k+1) * r_inc;
                        fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_reci_derivative_z_zz_cos(r,g);
                        res += r_inc*(fa+fb)/2.;
                }

                // 
        }
        // Bohr -> Angstrom
        return res*TO_BOHR_RADII;
}

////	Sine Part

//	Integral sx 'd/dgx'
double LonePairMatrix_H::reci_sx_grad_gx_sin( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double g )
{
        double res = 0.;
        double fa, fb, dr, mesh;
        double r,r_inc;
	// Unit 'g' Angs^-1 ... to Bohr^-1
	g = g*TO_BOHR_RADII;
        // Distance to Bohr

        for(int i=0;i<integral_knot.size()-1;i++)
        {
                dr   = integral_knot[i+1] - integral_knot[i];
                mesh = grid(dr);
                r_inc= dr/static_cast<double>(mesh);

		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
                for(int k=0;k<(int)mesh;k++)
                {
                        r  = integral_knot[i] + k * r_inc;
                        fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_reci_derivative_x_sx_sin(r,g);
                        r  = integral_knot[i] + (k+1) * r_inc;
                        fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_reci_derivative_x_sx_sin(r,g);
                        res += r_inc*(fa+fb)/2.;
                }

                // 
        }
        // Bohr -> Angstrom
        return res*TO_BOHR_RADII;
}

//	Integral sz 'd/dgz'
double LonePairMatrix_H::reci_sz_grad_gz_sin( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], double g )
{
        double res = 0.;
        double fa, fb, dr, mesh;
        double r,r_inc;
	// Unit 'g' Angs^-1 ... to Bohr^-1
	g = g*TO_BOHR_RADII;
        // Distance to Bohr

        for(int i=0;i<integral_knot.size()-1;i++)
        {
                dr   = integral_knot[i+1] - integral_knot[i];
                mesh = grid(dr);
                r_inc= dr/static_cast<double>(mesh);

		#pragma omp parallel for private(r,fa,fb) reduction(+:res) num_threads(NUM_OMP_THREAD) schedule(static)
                for(int k=0;k<(int)mesh;k++)
                {
                        r  = integral_knot[i] + k * r_inc;
                        fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_reci_derivative_z_sz_sin(r,g);
                        r  = integral_knot[i] + (k+1) * r_inc;
                        fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_reci_derivative_z_sz_sin(r,g);
                        res += r_inc*(fa+fb)/2.;
                }

                // 
        }
        // Bohr -> Angstrom
        return res*TO_BOHR_RADII;
}



//// DEV_TEST //// DEV_TEST
//// DEV_TEST //// DEV_TEST
//// DEV_TEST //// DEV_TEST
//// DEV_TEST //// DEV_TEST
//// DEV_TEST //// DEV_TEST
//// DEV_TEST //// DEV_TEST

void LonePairMatrix_H::test2()			// binary link test
{	std::cout << "LP mat Test 2\n";
}
double LonePairMatrix_H::NIntegral_test_real( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4] )
{
	double res, res_sz, res_xxyy, res_zz;
	res = res_sz = res_xxyy = res_zz = 0.;

	using std::cout;
	using std::endl;

	cout << "Knot stride :" << integral_knot.size() << endl;

	//Confirmed
	//for(int i=0;i<integral_knot.size();i++)
	//{     printf("%20.6e\n",integral_knot[i]);       }
	//for(int i=0;i<integral_knot.size()-1;i++)
	//{     printf("%20.6e%20.6e%20.6e%20.6e\n",Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i]);   }
	//for(int i=0;i<integral_knot.size()-1;i++)
	//{     printf("%20.6e%20.6e%20.6e%20.6e\n",Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i]);   }


	// private
	double as,bs,cs,ds;
	double ap,bp,cp,dp;
	double delta_r, incr_r;
	double fa,fb, r, ra, rb, rab;
	double fa_sz,fa_xxyy,fa_zz;
	double fb_sz,fb_xxyy,fb_zz;
	int grid;

	// shared
	//double w = 24.;				// approximately ~ 2048.
	//double w = 23484.;	// (1) - 0.5073390000 / (2) - 0.5635150000 / (3) - 0.6165529999
	double w = pow(2.,11);	// 10E-9 acc // (1) - 0.0524480001 
	double dist;
	double d;

	const double sig = 2.49475;
	double mm_e;
	double et;
	
/*
	r = integral_knot[0];
	for(int i=0;i<integral_knot.size()-1;i++)
	{	
		delta_r = integral_knot[i+1] - integral_knot[i];	// dr = r[i+1] - r[i];
		grid = ceil(fabs(w/log(delta_r)));
		//grid=2745;
		//if ( grid%2 != 0 ) { grid = grid+1; }
		//printf("%20.6e%20.6e\t%20.6d\n",delta_r,fabs(w/log(delta_r)),grid);
		printf("%20.12e\t%20.12e\t%20.12e\t%20.6d\n",r,delta_r,fabs(w/log(delta_r)),grid);
		r += delta_r;
	}
*/
	double pos_int = 0.;
	double k1,k2;
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		as = Rs[0][i];
		bs = Rs[1][i];
		cs = Rs[2][i];
		ds = Rs[3][i];

		ap = Rp[0][i];
		bp = Rp[1][i];
		cp = Rp[2][i];
		dp = Rp[3][i];

		k1 = integral_knot[i];
		k2 = integral_knot[i+1];

		// Position Integral of SLAM basis function ... confirmed Sep 09 2022 
		pos_int += (-0.5*(dp*ds*pow(k1,2)) - ((cs*dp + cp*ds)*pow(k1,3))/3. -
		((cp*cs + bs*dp + bp*ds)*pow(k1,4))/4. -
		((bs*cp + bp*cs + as*dp + ap*ds)*pow(k1,5))/5. -
		((bp*bs + as*cp + ap*cs)*pow(k1,6))/6. - ((as*bp + ap*bs)*pow(k1,7))/7. -
		(ap*as*pow(k1,8))/8. + (dp*ds*pow(k2,2))/2. + ((cs*dp + cp*ds)*pow(k2,3))/3. +
		((cp*cs + bs*dp + bp*ds)*pow(k2,4))/4. +
		((bs*cp + bp*cs + as*dp + ap*ds)*pow(k2,5))/5. +
		((bp*bs + as*cp + ap*cs)*pow(k2,6))/6. + ((as*bp + ap*bs)*pow(k2,7))/7. +
		(ap*as*pow(k2,8))/8.)/sqrt(3);
	}
	printf("PosInt Reference : %20.12lf\n",pos_int*TO_BOHR_RADII);
	cout << endl;

	et = -omp_get_wtime();
	//#pragma omp parallel private(as,ap,bs,bp,cs,cp,ds,dp,delta_r,incr_r,fa,fb,fa_sz,fb_sz,fa_xxyy,fb_xxyy,fa_zz,fb_zz,r,ra,rb,grid) num_threads(1)
	{
	
		dist = 0.05;
	for(int l=0;l<5000;l++)
	{
		dist = 0.005 + l*0.05;
		//dist = dist*1.0256;
		if( dist > 15. )
			break;
		//dist = 0.005 + l*0.005;
		d = dist;
		//mm_e = erfc(dist/sig)/dist;
	
		et = -omp_get_wtime();

#ifdef LP_DEBUG_1

		//#pragma omp for reduction(+:res,res_sz,res_xxyy,res_zz)
		for(int i=0;i<integral_knot.size()-1;i++)
		{
			// Below ... RealSpace NIntegrals
			delta_r = integral_knot[i+1] - integral_knot[i];	// dr = r[i+1] - r[i];
			grid = ceil(fabs(w/log(delta_r)));
			//if ( grid%2 != 0 ) { grid = grid+1; }
			//printf("%20.6e%20.6e\t%20.6d\n",delta_r,fabs(w/log(delta_r)),grid);

			as = Rs[0][i];
			bs = Rs[1][i];
			cs = Rs[2][i];
			ds = Rs[3][i];

			ap = Rp[0][i];
			bp = Rp[1][i];
			cp = Rp[2][i];
			dp = Rp[3][i];

			incr_r = delta_r/static_cast<double>(grid);

			for(int k=0;k<grid;k++)
			{
				ra = integral_knot[i] + k * incr_r;
				r  = ra;

				fa = pow(ds + cs*r + bs*pow(r,2.) + as*pow(r,3.),2.)*EnergyAngularIntegral_real_ss(sig,r,d);
				fa_sz   = (ds + cs*r + bs*pow(r,2.) + as*pow(r,3.))*(dp + cp*r + bp*pow(r,2.) + ap*pow(r,3.))*EnergyAngularIntegral_real_sz(sig,r,d);
				fa_xxyy = pow(dp + cp*r + bp*pow(r,2.) + ap*pow(r,3.),2.)*EnergyAngularIntegral_real_xx(sig,r,d);
				fa_zz   = pow(dp + cp*r + bp*pow(r,2.) + ap*pow(r,3.),2.)*EnergyAngularIntegral_real_zz(sig,r,d);

				rb = integral_knot[i] + (k+1) * incr_r;
				r  = rb;

				fb = pow(ds + cs*r + bs*pow(r,2.) + as*pow(r,3.),2.)*EnergyAngularIntegral_real_ss(sig,r,d);
				fb_sz   = (ds + cs*r + bs*pow(r,2.) + as*pow(r,3.))*(dp + cp*r + bp*pow(r,2.) + ap*pow(r,3.))*EnergyAngularIntegral_real_sz(sig,r,d);
				fb_xxyy = pow(dp + cp*r + bp*pow(r,2.) + ap*pow(r,3.),2.)*EnergyAngularIntegral_real_xx(sig,r,d);
				fb_zz   = pow(dp + cp*r + bp*pow(r,2.) + ap*pow(r,3.),2.)*EnergyAngularIntegral_real_zz(sig,r,d);

				//rab = rb - ra;
				res += incr_r*(fa+fb)/2.;
				res_sz   += incr_r*(fa_sz+fb_sz)/2.;
				res_xxyy += incr_r*(fa_xxyy+fb_xxyy)/2.;
				res_zz   += incr_r*(fa_zz+fb_zz)/2.;
				// Trapezoidal method
			}
		}
		et += omp_get_wtime();
		
		printf("%14.4lf\t%40.20e\t%40.20e\t%40.20e\t%40.20e\t%20.10lf\n",dist,res,res_sz,res_xxyy,res_zz,et);

		res = res_sz = res_xxyy = res_zz = 0;

#endif

#ifdef LP_DEBUG_2


/*
		res = res_xxyy = res_zz = 0.;

		for(int i=0;i<integral_knot.size()-1;i++)
		{
			//printf("%12.6lf%12.6lf\n",integral_knot[i],integral_knot[i+1]);
			res += EnergyIntegral_reci_ss(integral_knot[i],integral_knot[i+1],d,Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i]);
			res_xxyy += EnergyIntegral_reci_(integral_knot[i],integral_knot[i+1],d,Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i]);
			res_zz += EnergyIntegral_reci_zz(integral_knot[i],integral_knot[i+1],d,Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i]);

		}
		et += omp_get_wtime();

		printf("%14.4lf\t%40.20e\t%40.20e\t%40.20e\t%20.10lf\n",dist,res,res_xxyy,res_zz,et);
*/




		res = res_xxyy = res_zz = 0.;

		for(int i=0;i<integral_knot.size()-1;i++)
		{
			// Below ... RealSpace NIntegrals
			delta_r = integral_knot[i+1] - integral_knot[i];	// dr = r[i+1] - r[i];
			grid = ceil(fabs(w/log(delta_r)));
			//if ( grid%2 != 0 ) { grid = grid+1; }
			//printf("%20.6e%20.6e\t%20.6d\n",delta_r,fabs(w/log(delta_r)),grid);

			as = Rs[0][i];
			bs = Rs[1][i];
			cs = Rs[2][i];
			ds = Rs[3][i];

			ap = Rp[0][i];
			bp = Rp[1][i];
			cp = Rp[2][i];
			dp = Rp[3][i];

			incr_r = delta_r/static_cast<double>(grid);

			for(int k=0;k<grid;k++)
			{
				ra = integral_knot[i] + k * incr_r;
				r  = ra;

				fa = pow(ds + cs*r + bs*pow(r,2.) + as*pow(r,3.),2.)*EnergyAngularIntegral_reci_ss(r,d);
				fa_xxyy = pow(dp + cp*r + bp*pow(r,2.) + ap*pow(r,3.),2.)*EnergyAngularIntegral_reci_xx(r,d);
				fa_zz   = pow(dp + cp*r + bp*pow(r,2.) + ap*pow(r,3.),2.)*EnergyAngularIntegral_reci_zz(r,d);

				rb = integral_knot[i] + (k+1) * incr_r;
				r  = rb;

				fb = pow(ds + cs*r + bs*pow(r,2.) + as*pow(r,3.),2.)*EnergyAngularIntegral_reci_ss(r,d);//sin(d*r)/d/r;
				fb_xxyy = pow(dp + cp*r + bp*pow(r,2.) + ap*pow(r,3.),2.)*EnergyAngularIntegral_reci_xx(r,d);//3.*(-d*r*cos(d*r)+sin(d*r))/r/r/r/d/d/d;
				fb_zz   = pow(dp + cp*r + bp*pow(r,2.) + ap*pow(r,3.),2.)*EnergyAngularIntegral_reci_zz(r,d);//3.*(2.*d*r*cos(d*r)+(-2.+d*d*r*r)*sin(d*r))/d/d/d/r/r/r;

				res -= incr_r*(fa+fb)/2.;
				res_xxyy -= incr_r*(fa_xxyy+fb_xxyy)/2.;
				res_zz   -= incr_r*(fa_zz+fb_zz)/2.;
				// Trapezoidal method
			}
		}
		printf("up %14.4lf\t%40.20e\t%40.20e\t%40.20e\t%20.10lf\n",dist,res,res_xxyy,res_zz,et);

#endif
	}// l end

	}
	//et += omp_get_wtime();
	
	//printf("w: %20.6lf\t%40.12e\t%40.10lf\n",mm_e,res*HA_TO_EV_UNIT/2,et);

	return res;
}
