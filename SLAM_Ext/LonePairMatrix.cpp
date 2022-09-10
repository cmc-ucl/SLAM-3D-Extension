#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <omp.h>

#include "LonePairMatrix.hpp"

//#define LP_DEBUG_1

#define LP_DEBUG_2

// Internal Use Only
double radial( double a, double b, double c, double d, double r ) { return a*r*r*r + b*r*r + c*r + d; }
double grid( double delta_r ) { return ceil(fabs(GRID_W/log(delta_r))); }
//

const int LonePairMatrix::b_serach( const double dist, const std::vector<double>& integral_knot )
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

        if( (Rij(0) == 0) && (Rij(1) == 0) && (Rij(2) > 0) )                    // if vector 'Rij' is on z-axis
        {       this->transform_matrix(1,1) = 1.;
                this->transform_matrix(2,2) = 1.;
                this->transform_matrix(3,3) = 1.;                             // set lower-right block 3x3 matrix as I
        }
        else if( (Rij(0) == 0) && (Rij(1) == 0) && (Rij(2) < 0) )               // if vector 'Rij' is on negative z-axis
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

double LonePairMatrix_H::real_ss_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], const double sig, const double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	// Distance to Bohr
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		for(int k=0;k<mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_ss(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_ss(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res;
}

double LonePairMatrix_H::real_sz_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], const double sig, const double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	// Distance to Bohr
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		for(int k=0;k<mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_sz(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_sz(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res;
}

double LonePairMatrix_H::real_xx_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], const double sig, const double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	// Distance to Bohr
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		for(int k=0;k<mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_xx(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_xx(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res;
}

double LonePairMatrix_H::real_zz_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], const double sig, const double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	// Distance to Bohr
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		for(int k=0;k<mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_zz(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*radial(Rp[0][i],Rp[1][i],Rp[2][i],Rp[3][i],r)*EnergyAngularIntegral_real_zz(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res;
}

double LonePairMatrix_H::real_ss_grad_z_pc( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4], const double sig, const double d )
{
	double res = 0.;
	double fa, fb, dr, mesh;
	double r,r_inc;
	// Distance to Bohr
	
	for(int i=0;i<integral_knot.size()-1;i++)
	{
		dr   = integral_knot[i+1] - integral_knot[i];
		mesh = grid(dr);
		r_inc= dr/static_cast<double>(mesh);
		
		for(int k=0;k<mesh;k++)
		{
			r  = integral_knot[i] + k * r_inc;
			fa = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_derivative_z_ss(sig,r,d);
			r  = integral_knot[i] + (k+1) * r_inc;
			fb = radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*radial(Rs[0][i],Rs[1][i],Rs[2][i],Rs[3][i],r)*EnergyAngularIntegral_real_derivative_z_ss(sig,r,d);
			res += r_inc*(fa+fb)/2.;
		}
	}
	// eV Unit
	return res;
}


















void LonePairMatrix_H::test2()			// binary link test
{	std::cout << "LP mat Test 2\n";
}
const double LonePairMatrix_H::NIntegral_test_real( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4] )
{
	double res, res_sz, res_xxyy, res_zz;
	res = res_sz = res_xxyy = res_zz = 0.;

	using std::cout, std::endl;

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
