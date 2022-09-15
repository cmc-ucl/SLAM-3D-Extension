#ifndef __INTEGRAL_LIB
#define __INTEGRAL_LIB

double EnergyAngularIntegral_real_ss( double sig, double r, double d );
double EnergyAngularIntegral_real_sz( double sig, double r, double d );
double EnergyAngularIntegral_real_xx( double sig, double r, double d );		// xx = yy
double EnergyAngularIntegral_real_zz( double sig, double r, double d );

/* Note Sep 10 2022

 Derivative calculation is done for a function : f(|Ra + ra - Rb - Rb|) w.r.t. 'rb'

 i.e., df/d(rb;k) ... where k = x,y,z and the df/d(rb;k) result is equivalent with df/d(Rb;k)

 Consequently, df/d(rb;k) = df/d(Rb;k) = - df/d(Ra;k) ... means that can be used to calculate
 derivatives w.r.t. Ra interacting with a point charge, BUT the sign must be inversed !	*/
double EnergyAngularIntegral_real_derivative_x_sx( double sig, double r, double d );
double EnergyAngularIntegral_real_derivative_x_xz( double sig, double r, double d );

double EnergyAngularIntegral_real_derivative_z_ss( double sig, double r, double d );
double EnergyAngularIntegral_real_derivative_z_sz( double sig, double r, double d );
double EnergyAngularIntegral_real_derivative_z_xx( double sig, double r, double d );
double EnergyAngularIntegral_real_derivative_z_zz( double sig, double r, double d );

/* Note Sep 15 2022

 Derivative calculation is done in the same manner with the first derivatives

 i.e., d2f/d(rb;k)2 ... where k = x,y,z and the result is same with d2f/d(Rb;k)2 */
double EnergyAngularIntegral_real_derivative2_xx_ss( double sig, double r, double d );
double EnergyAngularIntegral_real_derivative2_xx_sz( double sig, double r, double d );
double EnergyAngularIntegral_real_derivative2_xx_xx( double sig, double r, double d );
double EnergyAngularIntegral_real_derivative2_xx_yy( double sig, double r, double d );
double EnergyAngularIntegral_real_derivative2_xx_zz( double sig, double r, double d );

double EnergyAngularIntegral_real_derivative2_xy_xy( double sig, double r, double d );

double EnergyAngularIntegral_real_derivative2_xz_sz( double sig, double r, double d );
double EnergyAngularIntegral_real_derivative2_xz_xz( double sig, double r, double d );

double EnergyAngularIntegral_real_derivative2_zz_ss( double sig, double r, double d );
double EnergyAngularIntegral_real_derivative2_zz_sz( double sig, double r, double d );
double EnergyAngularIntegral_real_derivative2_zz_xx( double sig, double r, double d );
double EnergyAngularIntegral_real_derivative2_zz_zz( double sig, double r, double d );

// Reciprocal Space Integrals
double EnergyAngularIntegral_reci_ss( double r, double g );
double EnergyAngularIntegral_reci_xx( double r, double g );
double EnergyAngularIntegral_reci_zz( double r, double g );
/*
double EnergyIntegral_reci_ss( double k1, double k2, double g, double as, double bs, double cs, double ds );
double EnergyIntegral_reci_xx( double k1, double k2, double g, double ap, double bp, double cp, double dp );
double EnergyIntegral_reci_zz( double k1, double k2, double g, double ap, double bp, double cp, double dp );
*/

#endif
