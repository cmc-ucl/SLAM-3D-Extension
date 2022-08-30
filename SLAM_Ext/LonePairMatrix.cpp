#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <omp.h>

#include "LonePairMatrix.hpp"

void LonePairMatrix_H::test2()			// binary link test
{	std::cout << "LP mat Test 2\n";
}

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
const double LonePairMatrix_H::NIntegral_test( const std::vector<double>& integral_knot, const std::vector<double> (&Rs)[4], const std::vector<double> (&Rp)[4] )
{
	double res=0;

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
	double delta_r, incr_r;
	double fa,fb, r, ra, rb, rab;
	double fa2,fb2;
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
	double res2=0;
	

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



	et = -omp_get_wtime();
	#pragma omp parallel private(as,bs,cs,ds,delta_r,incr_r,fa,fb,r,ra,rb,rab,grid) num_threads(1)
	{
	
	for(int l=0;l<1500;l++)
	{
		dist = 0.005 + l*0.005;
		d = dist;
		mm_e = erfc(dist/sig)/dist;
	
		et = -omp_get_wtime();
		#pragma omp for reduction(+:res,res2)
		for(int i=0;i<integral_knot.size()-1;i++)
		{
			delta_r = integral_knot[i+1] - integral_knot[i];	// dr = r[i+1] - r[i];
			grid = ceil(fabs(w/log(delta_r)));
			//grid=100;
			//if ( grid%2 != 0 ) { grid = grid+1; }
			//printf("%20.6e%20.6e\t%20.6d\n",delta_r,fabs(w/log(delta_r)),grid);

			as = Rs[0][i];
			bs = Rs[1][i];
			cs = Rs[2][i];
			ds = Rs[3][i];

			incr_r = delta_r/static_cast<double>(grid);

			for(int k=0;k<grid;k++)
			{
				ra = integral_knot[i] + k * incr_r;
				r  = ra;
				fa = (pow(ds + cs*r + bs*pow(r,2.) + as*pow(r,3.),2.)* 
				((2.*(-1. + pow(M_E,(4.*d*r)/pow(sig,2.)))*sig)/
				(pow(M_E,pow(d + r,2.)/pow(sig,2.))*sqrt(M_PI)) + 2.*(d + r)*erfc((d + r)/sig) - 
				2.*fabs(d - r)*erfc(fabs(d - r)/sig)))/(2.*d*r);
				//fa = pow(ds + cs*r + bs*pow(r,2.) + as*pow(r,3.),2.)*EnergyAngularIntegral_real_ss(sig,r,d);

				rb = integral_knot[i] + (k+1) * incr_r;
				r  = rb;
				//fb = pow(ds + cs*r + bs*pow(r,2.) + as*pow(r,3.),2.)*EnergyAngularIntegral_real_ss(sig,r,d);
				fb = (pow(ds + cs*r + bs*pow(r,2.) + as*pow(r,3.),2.)*
				((2.*(-1. + pow(M_E,(4.*d*r)/pow(sig,2.)))*sig)/
				(pow(M_E,pow(d + r,2.)/pow(sig,2.))*sqrt(M_PI)) + 2.*(d + r)*erfc((d + r)/sig) - 
				2.*fabs(d - r)*erfc(fabs(d - r)/sig)))/(2.*d*r);
				//rab = rb - ra;
				res += incr_r*(fa+fb)/2.;
				// Trapezoidal method
			}
		}
		et += omp_get_wtime();
		
	//std::cout << grid << std::endl;
	/*
		#pragma omp for reduction(+:res,res2)
		for(int i=0;i<integral_knot.size()-1;i++)
		{
			delta_r = integral_knot[i+1] - integral_knot[i];	// dr = r[i+1] - r[i];
			grid = ceil(fabs(w/log(delta_r)));
			//grid = 2746;
			if ( grid%2 != 0 ) { grid = grid+1; }
			//printf("%20.6e%20.6e\t%20.6d\n",delta_r,fabs(w/log(delta_r)),grid);

			as = Rs[0][i];
			bs = Rs[1][i];
			cs = Rs[2][i];
			ds = Rs[3][i];

			incr_r = delta_r/static_cast<double>(grid);

			for(int k=1;k<grid-1;k+=2)
			{
				r = integral_knot[i] + k * incr_r;

				res2+= 4.*incr_r/3.*
				(pow(ds + cs*r + bs*pow(r,2) + as*pow(r,3),2)*
				((2*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*sig)/
				(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)) + 2*(d + r)*erfc((d + r)/sig) - 
				2*fabs(d - r)*erfc(fabs(d - r)/sig)))/(2.*d*r);
			}
			for(int k=2;k<grid-2;k+=2)
			{
				r = integral_knot[i] + k * incr_r;

				res2+= 2.*incr_r/3.*
				(pow(ds + cs*r + bs*pow(r,2) + as*pow(r,3),2)*
				((2*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*sig)/
				(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)) + 2*(d + r)*erfc((d + r)/sig) - 
				2*fabs(d - r)*erfc(fabs(d - r)/sig)))/(2.*d*r);
			}

			r = integral_knot[i];
			res2 += incr_r/3.*
				(pow(ds + cs*r + bs*pow(r,2) + as*pow(r,3),2)*
				((2*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*sig)/
				(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)) + 2*(d + r)*erfc((d + r)/sig) - 
				2*fabs(d - r)*erfc(fabs(d - r)/sig)))/(2.*d*r);
			r = integral_knot[i+1];
			res2 += incr_r/3.*
				(pow(ds + cs*r + bs*pow(r,2) + as*pow(r,3),2)*
				((2*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*sig)/
				(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)) + 2*(d + r)*erfc((d + r)/sig) - 
				2*fabs(d - r)*erfc(fabs(d - r)/sig)))/(2.*d*r);

		}// i end	// knots
		et += omp_get_wtime();
	*/
		printf("d: %14.4lf\t%40.20e\t%40.20e\t%20.10lf\t%40.20e\t%40.20lf\n",dist,res,mm_e,et,res2,res-res2);

		res = 0;
		res2= 0;
	}// l end

	}
	//et += omp_get_wtime();
	
	//printf("w: %20.6lf\t%40.12e\t%40.10lf\n",mm_e,res*HA_TO_EV_UNIT/2,et);
	return res;
}
