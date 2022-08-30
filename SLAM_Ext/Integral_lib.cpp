#include <cmath>

double EnergyAngularIntegral_real_ss( double sig, double r, double d )
{
	return ((2*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*sig)/(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)) + 
	2*(d + r)*erfc((d + r)/sig) - 2*fabs(d - r)*erfc(fabs(d - r)/sig))/(4.*d*r);
}

double EnergyAngularIntegral_real_sz( double sig, double r, double d )
{
	return ((sig*(2*pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2))) + 
	2*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r + 
	(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(2*pow(r,2) - pow(sig,2))))/
	(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)) + 
	2*(pow(d,3) + pow(r,3))*erfc((d + r)/sig) - 
	3*(pow(d,2) + pow(r,2))*fabs(d - r)*erfc(fabs(d - r)/sig) + 
	pow(fabs(d - r),3)*erfc(fabs(d - r)/sig))/(4.*sqrt(3)*pow(d,2)*pow(r,2));
}

double EnergyAngularIntegral_real_xxyy( double sig, double r, double d )
{
	return (-3*((-2*sig*(pow(d,4)*(pow(M_E,pow(d - r,2)/pow(sig,2)) - 
	pow(M_E,pow(d + r,2)/pow(sig,2))) + 
	4*pow(d,3)*(pow(M_E,pow(d - r,2)/pow(sig,2)) + 
	pow(M_E,pow(d + r,2)/pow(sig,2)))*r + 
	4*d*(pow(M_E,pow(d - r,2)/pow(sig,2)) + pow(M_E,pow(d + r,2)/pow(sig,2)))*r*
	(pow(r,2) + pow(sig,2)) + 
	2*pow(d,2)*(pow(M_E,pow(d - r,2)/pow(sig,2)) - 
	pow(M_E,pow(d + r,2)/pow(sig,2)))*(3*pow(r,2) + pow(sig,2)) + 
	(pow(M_E,pow(d - r,2)/pow(sig,2)) - pow(M_E,pow(d + r,2)/pow(sig,2)))*
	(pow(r,4) + 2*pow(r,2)*pow(sig,2) + 2*pow(sig,4))))/
	(5.*pow(M_E,(2*(pow(d,2) + pow(r,2)))/pow(sig,2))*sqrt(M_PI)) + 
	(2*pow(d + r,5)*erfc((d + r)/sig))/5. - 
	(2*pow(fabs(d - r),5)*erfc(fabs(d - r)/sig))/5. + 
	pow(d,4)*((2*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*sig)/
	(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)) + 2*(d + r)*erfc((d + r)/sig) - 
	2*fabs(d - r)*erfc(fabs(d - r)/sig)) - 
	2*pow(d,2)*pow(r,2)*((2*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*sig)/
	(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)) + 2*(d + r)*erfc((d + r)/sig) - 
	2*fabs(d - r)*erfc(fabs(d - r)/sig)) + 
	pow(r,4)*((2*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*sig)/
	(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)) + 2*(d + r)*erfc((d + r)/sig) - 
	2*fabs(d - r)*erfc(fabs(d - r)/sig)) - 
	(4*pow(d,2)*((sig*(pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2))) - 
	2*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r + 
	(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(pow(r,2) + pow(sig,2))))/
	(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)) + 
	pow(d + r,3)*erfc((d + r)/sig) - pow(fabs(d - r),3)*erfc(fabs(d - r)/sig)))/3. - 
	(4*pow(r,2)*((sig*(pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2))) - 
	2*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r + 
	(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(pow(r,2) + pow(sig,2))))/
	(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)) + 
	pow(d + r,3)*erfc((d + r)/sig) - pow(fabs(d - r),3)*erfc(fabs(d - r)/sig)))/3.))/
	(32.*pow(d,3)*pow(r,3));
}

double EnergyAngularIntegral_real_zz( double sig, double r, double d )
{
	return (3*((-4*sig*(4*pow(d,4)*(pow(M_E,pow(d - r,2)/pow(sig,2)) - 
	pow(M_E,pow(d + r,2)/pow(sig,2))) - 
	4*pow(d,3)*(pow(M_E,pow(d - r,2)/pow(sig,2)) + 
	pow(M_E,pow(d + r,2)/pow(sig,2)))*r - 
	2*d*(pow(M_E,pow(d - r,2)/pow(sig,2)) + pow(M_E,pow(d + r,2)/pow(sig,2)))*r*
	(2*pow(r,2) - 3*pow(sig,2)) + 
	2*pow(d,2)*(pow(M_E,pow(d - r,2)/pow(sig,2)) - 
	pow(M_E,pow(d + r,2)/pow(sig,2)))*(7*pow(r,2) - pow(sig,2)) + 
	(pow(M_E,pow(d - r,2)/pow(sig,2)) - pow(M_E,pow(d + r,2)/pow(sig,2)))*
	(4*pow(r,4) - 2*pow(r,2)*pow(sig,2) + 3*pow(sig,4))))/
	(15.*pow(M_E,(2*(pow(d,2) + pow(r,2)))/pow(sig,2))*sqrt(M_PI)) + 
	(8*(2*pow(d,5) + 5*pow(d,3)*pow(r,2) + 5*pow(d,2)*pow(r,3) + 2*pow(r,5))*
	erfc((d + r)/sig))/15. - (8*(2*pow(d,4) + 2*pow(d,3)*r + 
	7*pow(d,2)*pow(r,2) + 2*d*pow(r,3) + 2*pow(r,4))*fabs(d - r)*
	erfc(fabs(d - r)/sig))/15.))/(16.*pow(d,3)*pow(r,3));
}
