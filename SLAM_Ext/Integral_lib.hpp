#ifndef __INTEGRAL_LIB
#define __INTEGRAL_LIB

double EnergyAngularIntegral_real_ss( double sig, double r, double d );
double EnergyAngularIntegral_real_sz( double sig, double r, double d );
double EnergyAngularIntegral_real_xxyy( double sig, double r, double d );
double EnergyAngularIntegral_real_zz( double sig, double r, double d );

double EnergyIntegral_reci_ss( double k1, double k2, double g, double as, double bs, double cs, double ds );
double EnergyIntegral_reci_xxyy( double k1, double k2, double g, double ap, double bp, double cp, double dp );
double EnergyIntegral_reci_zz( double k1, double k2, double g, double ap, double bp, double cp, double dp );

#endif
