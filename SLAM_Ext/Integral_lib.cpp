#include <cmath>
#include <gsl/gsl_sf_expint.h>

// Real Space Energy Integrals - Uses Numerical Trapzoidal Method
/*



*/

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

double EnergyAngularIntegral_real_xx( double sig, double r, double d )
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

////	////	////	////	////	////

// Energy Angular Integral 1st Derivatives

////	////	////	////	////	////

// derivative x - sx or y - sy ... usage : point charge derivative / lplp energy interaction

double EnergyAngularIntegral_real_derivative_x_sx( double sig, double r, double d )
{
	return ((sig*(2*pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2))) + 
	2*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r + 
	(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(2*pow(r,2) - pow(sig,2))))/
	(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)) + 
	2*(pow(d,3) + pow(r,3))*erfc((d + r)/sig) - 
	(2*pow(d - r,2)*(pow(d,2) + d*r + pow(r,2))*erfc(fabs(d - r)/sig))/fabs(d - r))/
	(4.*sqrt(3)*pow(d,3)*pow(r,2));
}

// derivative x - xz or y - yz ... usage : point charge derivative / lplp energy interaction

double EnergyAngularIntegral_real_derivative_x_xz( double sig, double r, double d )
{
	return (3*((sig*(4*pow(d,4)*(-1 + pow(M_E,(4*d*r)/pow(sig,2))) + 
	4*pow(d,3)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r + 
	2*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(2*pow(r,2) - 3*pow(sig,2)) + 
	2*pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(2*pow(r,2) - pow(sig,2)) + 
	(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*
	(4*pow(r,4) - 2*pow(r,2)*pow(sig,2) + 3*pow(sig,4))))/
	(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)) + 
	4*(pow(d,5) + pow(r,5))*erfc((d + r)/sig) - 
	(4*pow(d - r,2)*(pow(d,4) + pow(d,3)*r + pow(d,2)*pow(r,2) + d*pow(r,3) + 
	pow(r,4))*erfc(fabs(d - r)/sig))/fabs(d - r)))/(40.*pow(d,4)*pow(r,3));
}

// derivative z - ss           ... usage : point charge derivative / lplp energy interaction

double EnergyAngularIntegral_real_derivative_z_ss( double sig, double r, double d )
{
	return -0.5*(((-1 + pow(M_E,(4*d*r)/pow(sig,2)))*sig)/pow(M_E,pow(d + r,2)/pow(sig,2)) + 
	sqrt(M_PI)*r*erfc((d + r)/sig) + (sqrt(M_PI)*(d - r)*r*erfc(fabs(d - r)/sig))/fabs(d - r))/
	(pow(d,2)*sqrt(M_PI)*r);
}

// derivative z - sz           ... usage : point charge derivative / lplp energy interaction

double EnergyAngularIntegral_real_derivative_z_sz( double sig, double r, double d )
{
	return ((sig*(pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2))) - 
	2*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r - 
	(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(2*pow(r,2) - pow(sig,2))))/
	pow(M_E,pow(d + r,2)/pow(sig,2)) + 
	sqrt(M_PI)*(pow(d,3) - 2*pow(r,3))*erfc((d + r)/sig) - 
	(sqrt(M_PI)*(pow(d,4) - pow(d,3)*r + 2*d*pow(r,3) - 2*pow(r,4))*erfc(fabs(d - r)/sig))/
	fabs(d - r))/(2.*pow(d,3)*sqrt(3*M_PI)*pow(r,2));
}

// derivative z - xx = yy      ... usage : point charge derivative / lplp energy interaction

double EnergyAngularIntegral_real_derivative_z_xx( double sig, double r, double d )
{
	return ((sig*(-8*pow(d,4)*(-1 + pow(M_E,(4*d*r)/pow(sig,2))) - 
	8*pow(d,3)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r + 
	6*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(2*pow(r,2) - 3*pow(sig,2)) - 
	4*pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(2*pow(r,2) - pow(sig,2)) + 
	3*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*
	(4*pow(r,4) - 2*pow(r,2)*pow(sig,2) + 3*pow(sig,4))))/
	pow(M_E,pow(d + r,2)/pow(sig,2)) - 
	4*sqrt(M_PI)*pow(d + r,2)*(2*pow(d,3) - 4*pow(d,2)*r + 6*d*pow(r,2) - 3*pow(r,3))*
	erfc((d + r)/sig) + (4*sqrt(M_PI)*pow(d - r,3)*
	(2*pow(d,3) + 4*pow(d,2)*r + 6*d*pow(r,2) + 3*pow(r,3))*erfc(fabs(d - r)/sig))/
	fabs(d - r))/(40.*pow(d,4)*sqrt(M_PI)*pow(r,3));
}

// derivative z - zz           ... usage : point charge derivative / lplp energy interaction

double EnergyAngularIntegral_real_derivative_z_zz( double sig, double r, double d )
{
	return -0.05*((sig*(-8*pow(d,4)*(-1 + pow(M_E,(4*d*r)/pow(sig,2))) - 
	8*pow(d,3)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r + 
	6*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(2*pow(r,2) - 3*pow(sig,2)) + 
	2*pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(11*pow(r,2) + 2*pow(sig,2)) + 
	3*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*
	(4*pow(r,4) - 2*pow(r,2)*pow(sig,2) + 3*pow(sig,4))))/
	pow(M_E,pow(d + r,2)/pow(sig,2)) - 
	2*sqrt(M_PI)*(4*pow(d,5) - 5*pow(d,2)*pow(r,3) - 6*pow(r,5))*erfc((d + r)/sig) + 
	(2*sqrt(M_PI)*(4*pow(d,6) - 4*pow(d,5)*r + 5*pow(d,3)*pow(r,3) - 
	5*pow(d,2)*pow(r,4) + 6*d*pow(r,5) - 6*pow(r,6))*erfc(fabs(d - r)/sig))/
	fabs(d - r))/(pow(d,4)*sqrt(M_PI)*pow(r,3));
}

////	////	////	////	////	////

// Energy Angular Integral 2st Derivatives

////	////	////	////	////	////

double EnergyAngularIntegral_real_derivative2_xx_ss( double sig, double r, double d )
{
return ((-8*M_PI*(d*r + pow(r,2) - pow(d,2)*(-1 + sig))*pow(sig,3)*erfc((d + r)/sig))/(d + r) - 
(M_PI*pow(pow(d,2) - pow(r,2),2)*pow(sig,3)*erfc(fabs(d - r)/sig))/pow(fabs(d - r),3) + 
(2*M_PI*(3*pow(r,2) + pow(d,2)*(3 - 4*sig))*pow(sig,3)*erfc(fabs(d - r)/sig))/fabs(d - r) + 
3*M_PI*pow(sig,3)*fabs(d - r)*erfc(fabs(d - r)/sig) + sqrt(M_PI)*
(-((pow(sig,2)*(pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(-1 + sig) + 2*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(-1 + sig) + 
(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(pow(r,2)*(-1 + sig) + pow(sig,2)*(3 + 5*sig))))/pow(M_E,pow(d + r,2)/pow(sig,2))) - 
(-1 + sig)*(pow(d,4) + pow(r,4) + 6*pow(r,2)*pow(sig,2) - 2*pow(d,2)*(pow(r,2) - 3*pow(sig,2)))*
std::expint(-(pow(d - r,2)/pow(sig,2))) + (-1 + sig)*
(pow(d,4) + pow(r,4) + 6*pow(r,2)*pow(sig,2) - 2*pow(d,2)*(pow(r,2) - 3*pow(sig,2)))*
std::expint(-(pow(d + r,2)/pow(sig,2)))))/(16.*pow(d,3)*M_PI*r*pow(sig,4));
}
double EnergyAngularIntegral_real_derivative2_xx_sz( double sig, double r, double d )
{
return (sqrt(3)*((16*(-(d*pow(r,3)) - pow(r,4) + pow(d,4)*(-1 + sig) + pow(d,3)*r*(-1 + sig) + pow(d,2)*pow(r,2)*(-1 + sig))*pow(sig,3)*
erfc((d + r)/sig))/(d + r) - (pow(pow(d,2) - pow(r,2),2)*(pow(d,2) + pow(r,2))*pow(sig,3)*erfc(fabs(d - r)/sig))/
pow(fabs(d - r),3) - (pow(sig,3)*(-9*pow(r,4) + 2*pow(d,2)*pow(r,2)*(-3 + 4*sig) + pow(d,4)*(-9 + 8*sig))*
erfc(fabs(d - r)/sig))/fabs(d - r) + (9*pow(r,2) + pow(d,2)*(9 - 8*sig))*pow(sig,3)*fabs(d - r)*erfc(fabs(d - r)/sig) - 
pow(sig,3)*pow(fabs(d - r),3)*erfc(fabs(d - r)/sig) + 
((pow(sig,2)*(-(pow(d,4)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(-1 + sig)) - 
2*pow(d,3)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(-1 + sig) - 
2*pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(-1 + sig)*(pow(r,2) - 4*pow(sig,2)) - 
2*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(pow(r,2)*(-1 + sig) + pow(sig,2)*(1 + 7*sig)) - 
(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(pow(r,4)*(-1 + sig) + 8*pow(r,2)*pow(sig,2)*(1 + sig) - pow(sig,4)*(1 + 7*sig))))/
pow(M_E,pow(d + r,2)/pow(sig,2)) - (-1 + sig)*
(pow(d,6) + pow(r,6) + 9*pow(r,4)*pow(sig,2) - pow(d,4)*(pow(r,2) - 9*pow(sig,2)) - 
pow(d,2)*(pow(r,4) - 6*pow(r,2)*pow(sig,2)))*std::expint(-(pow(d - r,2)/pow(sig,2))) + 
(-1 + sig)*(pow(d,6) + pow(r,6) + 9*pow(r,4)*pow(sig,2) - pow(d,4)*(pow(r,2) - 9*pow(sig,2)) - 
pow(d,2)*(pow(r,4) - 6*pow(r,2)*pow(sig,2)))*std::expint(-(pow(d + r,2)/pow(sig,2))))/sqrt(M_PI)))/
(32.*pow(d,4)*pow(r,2)*pow(sig,4));
}
double EnergyAngularIntegral_real_derivative2_xx_xx( double sig, double r, double d )
{
return ((128*sqrt(M_PI)*pow(d - r,4)*pow(sig,3)*(-9*d*pow(r,3) - 9*pow(r,4) + pow(d,4)*(-9 + 5*sig) + pow(d,3)*r*(-9 + 5*sig) + pow(d,2)*pow(r,2)*(-9 + 5*sig))*erfc(fabs(d - r)/sig))/pow(fabs(d - r),3) + 
(pow(sig,2)*(45*pow(d,6)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(-1 + sig) + 90*pow(d,5)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(-1 + sig) + 
18*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(5*pow(r,4)*(-1 + sig) - 6*pow(sig,4)*(1 + 15*sig) + 2*pow(r,2)*pow(sig,2)*(7 + 25*sig)) + 
9*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(5*pow(r,6)*(-1 + sig) + 6*pow(sig,6)*(1 + 15*sig) - 2*pow(r,2)*pow(sig,4)*(7 + 25*sig) + pow(r,4)*pow(sig,2)*(73 + 55*sig)) - 
4*pow(d,3)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(45*pow(r,2)*(-1 + sig) - pow(sig,2)*(63 + 65*sig)) - 
pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(45*pow(r,4)*(-1 + sig) + 2*pow(sig,4)*(63 + 65*sig) - 2*pow(r,2)*pow(sig,2)*(171 + 85*sig)) - 
pow(d,4)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(45*pow(r,2)*(-1 + sig) + pow(sig,2)*(-657 + 145*sig)) - 
128*pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)*(d + r)*sig*(9*d*pow(r,3) - 9*pow(r,4) + pow(d,3)*r*(9 - 5*sig) + pow(d,4)*(-9 + 5*sig) + pow(d,2)*pow(r,2)*(-9 + 5*sig))*erfc((d + r)/sig)) + 
45*pow(M_E,pow(d + r,2)/pow(sig,2))*pow(pow(d,2) - pow(r,2),2)*(-1 + sig)*(pow(d,4) + pow(r,4) + 12*pow(r,2)*pow(sig,2) - 2*pow(d,2)*(pow(r,2) - 6*pow(sig,2)))*
(std::expint(-(pow(d - r,2)/pow(sig,2))) - std::expint(-(pow(d + r,2)/pow(sig,2)))))/pow(M_E,pow(d + r,2)/pow(sig,2)))/(1280.*pow(d,5)*sqrt(M_PI)*pow(r,3)*pow(sig,4));
}
double EnergyAngularIntegral_real_derivative2_xx_yy( double sig, double r, double d )
{
return ((15*pow(pow(d,2) - pow(r,2),4)*pow(sig,3)*erfc(fabs(d - r)/sig))/pow(fabs(d - r),3) + (60*pow(pow(d,2) - pow(r,2),2)*pow(sig,3)*(-3*pow(r,2) + pow(d,2)*(-3 + 4*sig))*erfc(fabs(d - r)/sig))/fabs(d - r)
 + 30*pow(sig,3)*(-9*pow(r,4) + 2*pow(d,2)*pow(r,2)*(-3 + 8*sig) + pow(d,4)*(-9 + 16*sig))*fabs(d - r)*erfc(fabs(d - r)/sig) - 
20*pow(sig,3)*(-3*pow(r,2) + pow(d,2)*(-3 + 4*sig))*pow(fabs(d - r),3)*erfc(fabs(d - r)/sig) - 9*pow(sig,3)*pow(fabs(d - r),5)*erfc(fabs(d - r)/sig) + 
(pow(sig,2)*(15*pow(d,6)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(-1 + sig) + 30*pow(d,5)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(-1 + sig) + 
6*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(5*pow(r,4)*(-1 + sig) - 6*pow(sig,4)*(1 + 15*sig) + 2*pow(r,2)*pow(sig,2)*(7 + 25*sig)) + 
3*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(5*pow(r,6)*(-1 + sig) + 6*pow(sig,6)*(1 + 15*sig) - 2*pow(r,2)*pow(sig,4)*(7 + 25*sig) + pow(r,4)*pow(sig,2)*(73 + 55*sig)) - 
4*pow(d,3)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(15*pow(r,2)*(-1 + sig) + pow(sig,2)*(-21 + 85*sig)) - 
pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(15*pow(r,4)*(-1 + sig) + 2*(21 - 85*sig)*pow(sig,4) + 2*pow(r,2)*pow(sig,2)*(-57 + 185*sig)) - 
pow(d,4)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(15*pow(r,2)*(-1 + sig) + pow(sig,2)*(-219 + 475*sig)) - 
128*pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)*(d + r)*sig*(3*d*pow(r,3) - 3*pow(r,4) + pow(d,3)*r*(3 - 5*sig) + pow(d,4)*(-3 + 5*sig) + pow(d,2)*pow(r,2)*(-3 + 5*sig))*erfc((d + r)/sig)) + 
15*pow(M_E,pow(d + r,2)/pow(sig,2))*pow(pow(d,2) - pow(r,2),2)*(-1 + sig)*(pow(d,4) + pow(r,4) + 12*pow(r,2)*pow(sig,2) - 2*pow(d,2)*(pow(r,2) - 6*pow(sig,2)))*
(std::expint(-(pow(d - r,2)/pow(sig,2))) - std::expint(-(pow(d + r,2)/pow(sig,2)))))/(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)))/(1280.*pow(d,5)*pow(r,3)*pow(sig,4));
}
double EnergyAngularIntegral_real_derivative2_xx_zz( double sig, double r, double d )
{
return ((32*pow(sig,3)*(-12*d*pow(r,5) - 12*pow(r,6) + 15*pow(d,4)*pow(r,2)*(-1 + sig) + 5*pow(d,3)*pow(r,3)*(-3 + 2*sig) + 5*pow(d,2)*pow(r,4)*(-3 + 2*sig) + 2*pow(d,6)*(-6 + 5*sig) + 
2*pow(d,5)*r*(-6 + 5*sig))*erfc((d + r)/sig))/(d + r) - (32*pow(d - r,2)*pow(sig,3)*(12*d*pow(r,5) - 12*pow(r,6) + 15*pow(d,4)*pow(r,2)*(-1 + sig) - 5*pow(d,3)*pow(r,3)*(-3 + 2*sig) + 
5*pow(d,2)*pow(r,4)*(-3 + 2*sig) + 2*pow(d,6)*(-6 + 5*sig) - 2*pow(d,5)*r*(-6 + 5*sig))*erfc(fabs(d - r)/sig))/pow(fabs(d - r),3) + 
(-((pow(sig,2)*(15*pow(d,6)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(-1 + sig) + 30*pow(d,5)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(-1 + sig) + 
pow(d,4)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(45*pow(r,2)*(-1 + sig) + (219 - 155*sig)*pow(sig,2)) + 4*pow(d,3)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(15*pow(r,2)*(-1 + sig) + (21 - 5*sig)*pow(sig,2)) + 
6*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(5*pow(r,4)*(-1 + sig) - 6*pow(sig,4)*(1 + 15*sig) + 2*pow(r,2)*pow(sig,2)*(7 + 25*sig)) + 
3*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(5*pow(r,6)*(-1 + sig) + 6*pow(sig,6)*(1 + 15*sig) - 2*pow(r,2)*pow(sig,4)*(7 + 25*sig) + pow(r,4)*pow(sig,2)*(73 + 55*sig)) + 
pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(45*pow(r,4)*(-1 + sig) + 2*pow(sig,4)*(-21 + 5*sig) + 2*pow(r,2)*pow(sig,2)*(147 + 125*sig))))/pow(M_E,pow(d + r,2)/pow(sig,2))) - 
15*(pow(d,2) + pow(r,2))*(-1 + sig)*(pow(d,6) - pow(d,2)*pow(r,4) + pow(r,6) + 12*pow(r,4)*pow(sig,2) - pow(d,4)*(pow(r,2) - 12*pow(sig,2)))*std::expint(-(pow(d - r,2)/pow(sig,2))) + 
15*(pow(d,2) + pow(r,2))*(-1 + sig)*(pow(d,6) - pow(d,2)*pow(r,4) + pow(r,6) + 12*pow(r,4)*pow(sig,2) - pow(d,4)*(pow(r,2) - 12*pow(sig,2)))*std::expint(-(pow(d + r,2)/pow(sig,2))))/sqrt(M_PI))/
(320.*pow(d,5)*pow(r,3)*pow(sig,4));
}
double EnergyAngularIntegral_real_derivative2_xy_xy( double sig, double r, double d )
{
return (3*((-128*sqrt(M_PI)*pow(d - r,4)*(pow(d,4) + pow(d,3)*r + pow(d,2)*pow(r,2) + d*pow(r,3) + pow(r,4))*pow(sig,3)*erfc(fabs(d - r)/sig))/pow(fabs(d - r),3) + 
(pow(sig,2)*(5*pow(d,6)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(-1 + sig) + 10*pow(d,5)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(-1 + sig) - 
4*pow(d,3)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(5*pow(r,2)*(-1 + sig) - pow(sig,2)*(7 + 25*sig)) + 
2*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(5*pow(r,4)*(-1 + sig) - 6*pow(sig,4)*(1 + 15*sig) + 2*pow(r,2)*pow(sig,2)*(7 + 25*sig)) - 
pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(5*pow(r,4)*(-1 + sig) + 2*pow(sig,4)*(7 + 25*sig) - 2*pow(r,2)*pow(sig,2)*(19 + 45*sig)) - 
pow(d,4)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(5*pow(r,2)*(-1 + sig) - pow(sig,2)*(73 + 55*sig)) + 
(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(5*pow(r,6)*(-1 + sig) + 6*pow(sig,6)*(1 + 15*sig) - 2*pow(r,2)*pow(sig,4)*(7 + 25*sig) + pow(r,4)*pow(sig,2)*(73 + 55*sig)) + 
128*pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)*(pow(d,5) + pow(r,5))*sig*erfc((d + r)/sig)) + 
5*pow(M_E,pow(d + r,2)/pow(sig,2))*pow(pow(d,2) - pow(r,2),2)*(-1 + sig)*(pow(d,4) + pow(r,4) + 12*pow(r,2)*pow(sig,2) - 2*pow(d,2)*(pow(r,2) - 6*pow(sig,2)))*
(std::expint(-(pow(d - r,2)/pow(sig,2))) - std::expint(-(pow(d + r,2)/pow(sig,2)))))/pow(M_E,pow(d + r,2)/pow(sig,2))))/(1280.*pow(d,5)*sqrt(M_PI)*pow(r,3)*pow(sig,4));
}
double EnergyAngularIntegral_real_derivative2_xz_sz( double sig, double r, double d )
{
return (3*((-16*pow(d - r,3)*pow(r,3)*pow(sig,3)*erfc(fabs(d - r)/sig))/pow(fabs(d - r),3) - (pow(sig,2)*
(-(pow(d,4)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(-1 + sig)) - 2*pow(d,3)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(-1 + sig) - 2*pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(-1 + sig)*pow(sig,2) + 
2*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(pow(r,2)*(-1 + sig) + pow(sig,2)*(1 + 7*sig)) + 
(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(pow(r,4)*(-1 + sig) + 8*pow(r,2)*pow(sig,2)*(1 + sig) - pow(sig,4)*(1 + 7*sig)) + 16*pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)*pow(r,3)*sig*erfc((d + r)/sig)) - 
pow(M_E,pow(d + r,2)/pow(sig,2))*(pow(d,2) - pow(r,2))*(-1 + sig)*(pow(d,4) + pow(r,4) + 9*pow(r,2)*pow(sig,2) + pow(d,2)*(-2*pow(r,2) + 3*pow(sig,2)))*
(std::expint(-(pow(d - r,2)/pow(sig,2))) - std::expint(-(pow(d + r,2)/pow(sig,2)))))/(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI))))/(32.*pow(d,4)*pow(r,2)*pow(sig,4));
}
double EnergyAngularIntegral_real_derivative2_xz_xz( double sig, double r, double d )
{
return (3*((5*pow(pow(d,2) - pow(r,2),3)*(pow(d,2) + pow(r,2))*pow(sig,3)*erfc(fabs(d - r)/sig))/pow(fabs(d - r),3) - (30*(pow(d,6) + pow(d,2)*pow(r,4) - 
2*pow(r,6))*pow(sig,3)*erfc(fabs(d - r)/sig))/fabs(d - r) + 30*pow(r,2)*(pow(d,2) + 3*pow(r,2))*pow(sig,3)*fabs(d - r)*erfc(fabs(d - r)/sig) - 10*(pow(d,2) + 
2*pow(r,2))*pow(sig,3)*pow(fabs(d - r),3)*erfc(fabs(d - r)/sig) + 3*pow(sig,3)*pow(fabs(d - r),5)*erfc(fabs(d - r)/sig) + (pow(sig,2)*(5*pow(d,6)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(-1 + sig) + 
10*pow(d,5)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(-1 + sig) + 8*pow(d,3)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*pow(sig,2)*(-1 + 5*sig) + pow(d,4)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(5*pow(r,2)*(-1 + sig) + 
pow(sig,2)*(7 + 25*sig)) - 2*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(5*pow(r,4)*(-1 + sig) - 6*pow(sig,4)*(1 + 15*sig) + 2*pow(r,2)*pow(sig,2)*(7 + 25*sig)) - 
pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(5*pow(r,4)*(-1 + sig) + 4*pow(sig,4)*(-1 + 5*sig) + 2*pow(r,2)*pow(sig,2)*(9 + 55*sig)) - 
(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(5*pow(r,6)*(-1 + sig) + 6*pow(sig,6)*(1 + 15*sig) - 2*pow(r,2)*pow(sig,4)*(7 + 25*sig) + pow(r,4)*pow(sig,2)*(73 + 55*sig)) + 
32*pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)*(pow(d,5) - 4*pow(r,5))*sig*erfc((d + r)/sig)) + 
5*pow(M_E,pow(d + r,2)/pow(sig,2))*(pow(d,2) - pow(r,2))*(-1 + sig)*(pow(d,6) + pow(r,6) + 12*pow(r,4)*pow(sig,2) - pow(d,4)*(pow(r,2) - 6*pow(sig,2)) - 
pow(d,2)*(pow(r,4) - 6*pow(r,2)*pow(sig,2)))*(std::expint(-(pow(d - r,2)/pow(sig,2))) - std::expint(-(pow(d + r,2)/pow(sig,2)))))/(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI))))/
(320.*pow(d,5)*pow(r,3)*pow(sig,4));
}
double EnergyAngularIntegral_real_derivative2_zz_ss( double sig, double r, double d )
{
return -(-((pow(sig,2)*(2*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(-1 + sig) + pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(-1 + 9*sig) + (-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(pow(r,2)*(-1 + sig) + 
pow(sig,2)*(3 + 5*sig))))/pow(M_E,pow(d + r,2)/pow(sig,2))) - (4*sqrt(M_PI)*(2*d*r + 2*pow(r,2) + pow(d,2)*(-1 + sig))*pow(sig,3)*erfc((d + r)/sig))/(d + r) - 
(sqrt(M_PI)*pow(pow(d,2) - pow(r,2),2)*pow(sig,3)*erfc(fabs(d - r)/sig))/pow(fabs(d - r),3) + (2*sqrt(M_PI)*pow(sig,3)*(3*pow(r,2) + pow(d,2)*(-3 + 2*sig))*erfc(fabs(d - r)/sig))/fabs(d - r) + 
3*sqrt(M_PI)*pow(sig,3)*fabs(d - r)*erfc(fabs(d - r)/sig) - (pow(d,2) - pow(r,2))*(-1 + sig)*(pow(d,2) - pow(r,2) - 6*pow(sig,2))*std::expint(-(pow(d - r,2)/pow(sig,2))) + 
(pow(d,2) - pow(r,2))*(-1 + sig)*(pow(d,2) - pow(r,2) - 6*pow(sig,2))*std::expint(-(pow(d + r,2)/pow(sig,2))))/(8.*pow(d,3)*sqrt(M_PI)*r*pow(sig,4));
}
double EnergyAngularIntegral_real_derivative2_zz_sz( double sig, double r, double d )
{
return (sqrt(3)*((8*(2*d*pow(r,3) + 2*pow(r,4) + pow(d,4)*(-1 + sig) + pow(d,3)*r*(-1 + sig) + pow(d,2)*pow(r,2)*(-1 + sig))*pow(sig,3)*erfc((d + r)/sig))/(d + r) + 
(pow(pow(d,2) - pow(r,2),2)*(pow(d,2) + pow(r,2))*pow(sig,3)*erfc(fabs(d - r)/sig))/pow(fabs(d - r),3) + 
((-9*pow(r,4) + pow(d,4)*(3 - 4*sig) + 2*pow(d,2)*pow(r,2)*(3 - 2*sig))*pow(sig,3)*erfc(fabs(d - r)/sig))/fabs(d - r) - pow(sig,3)*(9*pow(r,2) + pow(d,2)*(-3 + 4*sig))*fabs(d - r)*erfc(fabs(d - r)/sig) + 
pow(sig,3)*pow(fabs(d - r),3)*erfc(fabs(d - r)/sig) + ((pow(sig,2)*(pow(d,4)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(-1 + sig) + 2*pow(d,3)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(-1 + 9*sig) + 
2*pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(pow(r,2)*(-1 + sig) - 2*pow(sig,2)*(1 + sig)) + 2*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(pow(r,2)*(-1 + sig) + pow(sig,2)*(1 + 7*sig)) + 
(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(pow(r,4)*(-1 + sig) + 8*pow(r,2)*pow(sig,2)*(1 + sig) - pow(sig,4)*(1 + 7*sig))))/pow(M_E,pow(d + r,2)/pow(sig,2)) + 
(pow(d,2) - pow(r,2))*(-1 + sig)*(pow(d,4) - 3*pow(d,2)*pow(sig,2) - pow(r,2)*(pow(r,2) + 9*pow(sig,2)))*std::expint(-(pow(d - r,2)/pow(sig,2))) + 
(pow(d,2) - pow(r,2))*(-1 + sig)*(-pow(d,4) + pow(r,4) + 3*pow(d,2)*pow(sig,2) + 9*pow(r,2)*pow(sig,2))*std::expint(-(pow(d + r,2)/pow(sig,2))))/sqrt(M_PI)))/(16.*pow(d,4)*pow(r,2)*pow(sig,4));
}
double EnergyAngularIntegral_real_derivative2_zz_xx( double sig, double r, double d )
{
return ((-15*pow(pow(d,2) - pow(r,2),4)*pow(sig,3)*erfc(fabs(d - r)/sig))/pow(fabs(d - r),3) + (60*pow(pow(d,2) - pow(r,2),2)*pow(sig,3)*(3*pow(r,2) + pow(d,2)*sig)*erfc(fabs(d - r)/sig))/fabs(d - r) + 
30*pow(sig,3)*(9*pow(r,4) + 2*pow(d,2)*pow(r,2)*(-3 + 2*sig) + pow(d,4)*(-3 + 4*sig))*fabs(d - r)*erfc(fabs(d - r)/sig) - 20*pow(sig,3)*(3*pow(r,2) + pow(d,2)*sig)*pow(fabs(d - r),3)*erfc(fabs(d - r)/sig) + 
9*pow(sig,3)*pow(fabs(d - r),5)*erfc(fabs(d - r)/sig) + (pow(sig,2)*(-15*pow(d,6)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(-1 + sig) - 30*pow(d,5)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(-1 + sig) + 
pow(d,4)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(15*pow(r,2)*(-1 + sig) + (81 - 145*sig)*pow(sig,2)) - 
6*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(5*pow(r,4)*(-1 + sig) - 6*pow(sig,4)*(1 + 15*sig) + 2*pow(r,2)*pow(sig,2)*(7 + 25*sig)) - 
3*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(5*pow(r,6)*(-1 + sig) + 6*pow(sig,6)*(1 + 15*sig) - 2*pow(r,2)*pow(sig,4)*(7 + 25*sig) + pow(r,4)*pow(sig,2)*(73 + 55*sig)) + 
4*pow(d,3)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(15*pow(r,2)*(-1 + sig) + pow(sig,2)*(9 + 95*sig)) + 
pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(15*pow(r,4)*(-1 + sig) + 2*pow(r,2)*(93 - 125*sig)*pow(sig,2) - 2*pow(sig,4)*(9 + 95*sig)) - 
32*pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)*(d + r)*sig*(-12*d*pow(r,3) + 12*pow(r,4) + pow(d,3)*r*(3 - 5*sig) + pow(d,4)*(-3 + 5*sig) + pow(d,2)*pow(r,2)*(-3 + 5*sig))*erfc((d + r)/sig) - 
15*pow(M_E,pow(d + r,2)/pow(sig,2))*pow(pow(d,2) - pow(r,2),2)*(-1 + sig)*(pow(d,4) - 2*pow(d,2)*pow(r,2) + pow(r,4) + 12*pow(r,2)*pow(sig,2))*
(std::expint(-(pow(d - r,2)/pow(sig,2))) - std::expint(-(pow(d + r,2)/pow(sig,2))))))/(pow(M_E,pow(d + r,2)/pow(sig,2))*sqrt(M_PI)))/(320.*pow(d,5)*pow(r,3)*pow(sig,4));
}
double EnergyAngularIntegral_real_derivative2_zz_zz( double sig, double r, double d )
{
return ((16*pow(sig,3)*(24*d*pow(r,5) + 24*pow(r,6) + 15*pow(d,4)*pow(r,2)*(-1 + sig) + 10*pow(d,3)*pow(r,3)*sig + 10*pow(d,2)*pow(r,4)*sig + 2*pow(d,6)*(-3 + 5*sig) + 
2*pow(d,5)*r*(-3 + 5*sig))*erfc((d + r)/sig))/(d + r) + (15*pow(pow(d,4) - pow(r,4),2)*pow(sig,3)*erfc(fabs(d - r)/sig))/pow(fabs(d - r),3) - 
(60*(pow(d,2) + pow(r,2))*pow(sig,3)*(3*pow(r,4) + pow(d,2)*pow(r,2)*(-3 + sig) + pow(d,4)*sig)*erfc(fabs(d - r)/sig))/fabs(d - r) - 
30*pow(sig,3)*(9*pow(r,4) + 4*pow(d,2)*pow(r,2)*sig + pow(d,4)*(-3 + 4*sig))*fabs(d - r)*erfc(fabs(d - r)/sig) + 20*pow(sig,3)*(3*pow(r,2) + pow(d,2)*sig)*pow(fabs(d - r),3)*erfc(fabs(d - r)/sig) - 
9*pow(sig,3)*pow(fabs(d - r),5)*erfc(fabs(d - r)/sig) + ((pow(sig,2)*(15*pow(d,6)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(-1 + sig) + 30*pow(d,5)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(-1 + sig) + 
6*d*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(5*pow(r,4)*(-1 + sig) - 6*pow(sig,4)*(1 + 15*sig) + 2*pow(r,2)*pow(sig,2)*(7 + 25*sig)) + 
3*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(5*pow(r,6)*(-1 + sig) + 6*pow(sig,6)*(1 + 15*sig) - 2*pow(r,2)*pow(sig,4)*(7 + 25*sig) + pow(r,4)*pow(sig,2)*(73 + 55*sig)) + 
4*pow(d,3)*(1 + pow(M_E,(4*d*r)/pow(sig,2)))*r*(15*pow(r,2)*(-1 + sig) - pow(sig,2)*(9 + 95*sig)) + pow(d,4)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(15*pow(r,2)*(-3 + 35*sig) + pow(sig,2)*(-81 + 145*sig)) + 
pow(d,2)*(-1 + pow(M_E,(4*d*r)/pow(sig,2)))*(45*pow(r,4)*(-1 + sig) + 2*pow(sig,4)*(9 + 95*sig) + 2*pow(r,2)*pow(sig,2)*(-3 + 275*sig))))/pow(M_E,pow(d + r,2)/pow(sig,2)) + 
15*(pow(d,4) - pow(r,4))*(-1 + sig)*(pow(d,4) - pow(r,2)*(pow(r,2) + 12*pow(sig,2)))*std::expint(-(pow(d - r,2)/pow(sig,2))) - 
15*(pow(d,4) - pow(r,4))*(-1 + sig)*(pow(d,4) - pow(r,2)*(pow(r,2) + 12*pow(sig,2)))*std::expint(-(pow(d + r,2)/pow(sig,2))))/sqrt(M_PI))/(160.*pow(d,5)*pow(r,3)*pow(sig,4));
}



////	////	////	////	////	////

// Energy Reciprocal space integrals

////	////	////	////	////	////

double EnergyAngularIntegral_reci_ss( const double r, const double g )
{
return sin(g*r)/g/r;
}

double EnergyAngularIntegral_reci_xx( const double r, const double g )
{
return 3.*(-g*r*cos(g*r)+sin(g*r))/g/g/g/r/r/r;
}

double EnergyAngularIntegral_reci_zz( const double r, const double g )
{
return 3.*(2.*g*r*cos(g*r)+(-2.+g*g*r*r)*sin(g*r))/g/g/g/r/r/r;
}

/*

double EnergyIntegral_reci_ss( double k1, double k2, double g, double as, double bs, double cs, double ds )
{
return (48*as*bs*g*cos(g*k1) - 4*bs*cs*pow(g,3)*cos(g*k1) - 4*as*ds*pow(g,3)*cos(g*k1) + 
2*cs*ds*pow(g,5)*cos(g*k1) + 120*pow(as,2)*g*k1*cos(g*k1) - 6*pow(bs,2)*pow(g,3)*k1*cos(g*k1) - 
12*as*cs*pow(g,3)*k1*cos(g*k1) + pow(cs,2)*pow(g,5)*k1*cos(g*k1) + 
2*bs*ds*pow(g,5)*k1*cos(g*k1) - 24*as*bs*pow(g,3)*pow(k1,2)*cos(g*k1) + 
2*bs*cs*pow(g,5)*pow(k1,2)*cos(g*k1) + 2*as*ds*pow(g,5)*pow(k1,2)*cos(g*k1) - 
20*pow(as,2)*pow(g,3)*pow(k1,3)*cos(g*k1) + pow(bs,2)*pow(g,5)*pow(k1,3)*cos(g*k1) + 
2*as*cs*pow(g,5)*pow(k1,3)*cos(g*k1) + 2*as*bs*pow(g,5)*pow(k1,4)*cos(g*k1) + 
pow(as,2)*pow(g,5)*pow(k1,5)*cos(g*k1) - 48*as*bs*g*cos(g*k2) + 4*bs*cs*pow(g,3)*cos(g*k2) + 
4*as*ds*pow(g,3)*cos(g*k2) - 2*cs*ds*pow(g,5)*cos(g*k2) - 120*pow(as,2)*g*k2*cos(g*k2) + 
6*pow(bs,2)*pow(g,3)*k2*cos(g*k2) + 12*as*cs*pow(g,3)*k2*cos(g*k2) - 
pow(cs,2)*pow(g,5)*k2*cos(g*k2) - 2*bs*ds*pow(g,5)*k2*cos(g*k2) + 
24*as*bs*pow(g,3)*pow(k2,2)*cos(g*k2) - 2*bs*cs*pow(g,5)*pow(k2,2)*cos(g*k2) - 
2*as*ds*pow(g,5)*pow(k2,2)*cos(g*k2) + 20*pow(as,2)*pow(g,3)*pow(k2,3)*cos(g*k2) - 
pow(bs,2)*pow(g,5)*pow(k2,3)*cos(g*k2) - 2*as*cs*pow(g,5)*pow(k2,3)*cos(g*k2) - 
2*as*bs*pow(g,5)*pow(k2,4)*cos(g*k2) - pow(as,2)*pow(g,5)*pow(k2,5)*cos(g*k2) - 
120*pow(as,2)*sin(g*k1) + 6*pow(bs,2)*pow(g,2)*sin(g*k1) + 12*as*cs*pow(g,2)*sin(g*k1) - 
pow(cs,2)*pow(g,4)*sin(g*k1) - 2*bs*ds*pow(g,4)*sin(g*k1) + 48*as*bs*pow(g,2)*k1*sin(g*k1) - 
4*bs*cs*pow(g,4)*k1*sin(g*k1) - 4*as*ds*pow(g,4)*k1*sin(g*k1) + 
60*pow(as,2)*pow(g,2)*pow(k1,2)*sin(g*k1) - 3*pow(bs,2)*pow(g,4)*pow(k1,2)*sin(g*k1) - 
6*as*cs*pow(g,4)*pow(k1,2)*sin(g*k1) - 8*as*bs*pow(g,4)*pow(k1,3)*sin(g*k1) - 
5*pow(as,2)*pow(g,4)*pow(k1,4)*sin(g*k1) + 120*pow(as,2)*sin(g*k2) - 
6*pow(bs,2)*pow(g,2)*sin(g*k2) - 12*as*cs*pow(g,2)*sin(g*k2) + pow(cs,2)*pow(g,4)*sin(g*k2) + 
2*bs*ds*pow(g,4)*sin(g*k2) - 48*as*bs*pow(g,2)*k2*sin(g*k2) + 4*bs*cs*pow(g,4)*k2*sin(g*k2) + 
4*as*ds*pow(g,4)*k2*sin(g*k2) - 60*pow(as,2)*pow(g,2)*pow(k2,2)*sin(g*k2) + 
3*pow(bs,2)*pow(g,4)*pow(k2,2)*sin(g*k2) + 6*as*cs*pow(g,4)*pow(k2,2)*sin(g*k2) + 
8*as*bs*pow(g,4)*pow(k2,3)*sin(g*k2) + 5*pow(as,2)*pow(g,4)*pow(k2,4)*sin(g*k2) - 
pow(ds,2)*pow(g,6)*gsl_sf_Si(g*k1) + pow(ds,2)*pow(g,6)*gsl_sf_Si(g*k2))/pow(g,7);
}

double EnergyIntegral_reci_xxyy( double k1, double k2, double g, double ap, double bp, double cp, double dp )
{
return (3*((-16*ap*bp*cos(g*k1))/pow(g,3) + (4*bp*cp*cos(g*k1))/g + (4*ap*dp*cos(g*k1))/g - 
(pow(dp,2)*g*cos(g*k1))/(2.*k1) - (30*pow(ap,2)*k1*cos(g*k1))/pow(g,3) + 
(3*pow(bp,2)*k1*cos(g*k1))/g + (6*ap*cp*k1*cos(g*k1))/g + (8*ap*bp*pow(k1,2)*cos(g*k1))/g + 
(5*pow(ap,2)*pow(k1,3)*cos(g*k1))/g + (16*ap*bp*cos(g*k2))/pow(g,3) - (4*bp*cp*cos(g*k2))/g - 
	(4*ap*dp*cos(g*k2))/g + (pow(dp,2)*g*cos(g*k2))/(2.*k2) + 
	(30*pow(ap,2)*k2*cos(g*k2))/pow(g,3) - (3*pow(bp,2)*k2*cos(g*k2))/g - 
	(6*ap*cp*k2*cos(g*k2))/g - (8*ap*bp*pow(k2,2)*cos(g*k2))/g - 
	(5*pow(ap,2)*pow(k2,3)*cos(g*k2))/g + pow(cp,2)*sin(g*k1) + 2*bp*dp*sin(g*k1) + 
	(30*pow(ap,2)*sin(g*k1))/pow(g,4) - (3*pow(bp,2)*sin(g*k1))/pow(g,2) - 
	(6*ap*cp*sin(g*k1))/pow(g,2) + (pow(dp,2)*sin(g*k1))/(2.*pow(k1,2)) + (2*cp*dp*sin(g*k1))/k1 + 
	2*bp*cp*k1*sin(g*k1) + 2*ap*dp*k1*sin(g*k1) - (16*ap*bp*k1*sin(g*k1))/pow(g,2) + 
	pow(bp,2)*pow(k1,2)*sin(g*k1) + 2*ap*cp*pow(k1,2)*sin(g*k1) - 
	(15*pow(ap,2)*pow(k1,2)*sin(g*k1))/pow(g,2) + 2*ap*bp*pow(k1,3)*sin(g*k1) + 
	pow(ap,2)*pow(k1,4)*sin(g*k1) - pow(cp,2)*sin(g*k2) - 2*bp*dp*sin(g*k2) - 
	(30*pow(ap,2)*sin(g*k2))/pow(g,4) + (3*pow(bp,2)*sin(g*k2))/pow(g,2) + 
	(6*ap*cp*sin(g*k2))/pow(g,2) - (pow(dp,2)*sin(g*k2))/(2.*pow(k2,2)) - (2*cp*dp*sin(g*k2))/k2 - 
	2*bp*cp*k2*sin(g*k2) - 2*ap*dp*k2*sin(g*k2) + (16*ap*bp*k2*sin(g*k2))/pow(g,2) - 
	pow(bp,2)*pow(k2,2)*sin(g*k2) - 2*ap*cp*pow(k2,2)*sin(g*k2) + 
	(15*pow(ap,2)*pow(k2,2)*sin(g*k2))/pow(g,2) - 2*ap*bp*pow(k2,3)*sin(g*k2) - 
	pow(ap,2)*pow(k2,4)*sin(g*k2) - ((2*pow(cp,2) + dp*(4*bp + dp*pow(g,2)))*gsl_sf_Si(g*k1))/2. + 
	(pow(cp,2) + 2*bp*dp + (pow(dp,2)*pow(g,2))/2.)*gsl_sf_Si(g*k2)))/pow(g,3);
}


double EnergyIntegral_reci_zz( double k1, double k2, double g, double ap, double bp, double cp, double dp )
{
	return (3*((80*ap*bp*cos(g*k1))/pow(g,3) - (12*bp*cp*cos(g*k1))/g - (12*ap*dp*cos(g*k1))/g + 
	2*cp*dp*g*cos(g*k1) + (pow(dp,2)*g*cos(g*k1))/k1 + (180*pow(ap,2)*k1*cos(g*k1))/pow(g,3) - 
	(12*pow(bp,2)*k1*cos(g*k1))/g - (24*ap*cp*k1*cos(g*k1))/g + pow(cp,2)*g*k1*cos(g*k1) + 
	2*bp*dp*g*k1*cos(g*k1) - (40*ap*bp*pow(k1,2)*cos(g*k1))/g + 2*bp*cp*g*pow(k1,2)*cos(g*k1) + 
	2*ap*dp*g*pow(k1,2)*cos(g*k1) - (30*pow(ap,2)*pow(k1,3)*cos(g*k1))/g + 
	pow(bp,2)*g*pow(k1,3)*cos(g*k1) + 2*ap*cp*g*pow(k1,3)*cos(g*k1) + 
	2*ap*bp*g*pow(k1,4)*cos(g*k1) + pow(ap,2)*g*pow(k1,5)*cos(g*k1) - 
	(80*ap*bp*cos(g*k2))/pow(g,3) + (12*bp*cp*cos(g*k2))/g + (12*ap*dp*cos(g*k2))/g - 
	2*cp*dp*g*cos(g*k2) - (pow(dp,2)*g*cos(g*k2))/k2 - (180*pow(ap,2)*k2*cos(g*k2))/pow(g,3) + 
	(12*pow(bp,2)*k2*cos(g*k2))/g + (24*ap*cp*k2*cos(g*k2))/g - pow(cp,2)*g*k2*cos(g*k2) - 
	2*bp*dp*g*k2*cos(g*k2) + (40*ap*bp*pow(k2,2)*cos(g*k2))/g - 2*bp*cp*g*pow(k2,2)*cos(g*k2) - 
	2*ap*dp*g*pow(k2,2)*cos(g*k2) + (30*pow(ap,2)*pow(k2,3)*cos(g*k2))/g - 
	pow(bp,2)*g*pow(k2,3)*cos(g*k2) - 2*ap*cp*g*pow(k2,3)*cos(g*k2) - 
	2*ap*bp*g*pow(k2,4)*cos(g*k2) - pow(ap,2)*g*pow(k2,5)*cos(g*k2) - 3*pow(cp,2)*sin(g*k1) - 
	6*bp*dp*sin(g*k1) - (180*pow(ap,2)*sin(g*k1))/pow(g,4) + (12*pow(bp,2)*sin(g*k1))/pow(g,2) + 
	(24*ap*cp*sin(g*k1))/pow(g,2) - (pow(dp,2)*sin(g*k1))/pow(k1,2) - (4*cp*dp*sin(g*k1))/k1 - 
	8*bp*cp*k1*sin(g*k1) - 8*ap*dp*k1*sin(g*k1) + (80*ap*bp*k1*sin(g*k1))/pow(g,2) - 
	5*pow(bp,2)*pow(k1,2)*sin(g*k1) - 10*ap*cp*pow(k1,2)*sin(g*k1) + 
	(90*pow(ap,2)*pow(k1,2)*sin(g*k1))/pow(g,2) - 12*ap*bp*pow(k1,3)*sin(g*k1) - 
	7*pow(ap,2)*pow(k1,4)*sin(g*k1) + 3*pow(cp,2)*sin(g*k2) + 6*bp*dp*sin(g*k2) + 
	(180*pow(ap,2)*sin(g*k2))/pow(g,4) - (12*pow(bp,2)*sin(g*k2))/pow(g,2) - 
	(24*ap*cp*sin(g*k2))/pow(g,2) + (pow(dp,2)*sin(g*k2))/pow(k2,2) + (4*cp*dp*sin(g*k2))/k2 + 
	8*bp*cp*k2*sin(g*k2) + 8*ap*dp*k2*sin(g*k2) - (80*ap*bp*k2*sin(g*k2))/pow(g,2) + 
	5*pow(bp,2)*pow(k2,2)*sin(g*k2) + 10*ap*cp*pow(k2,2)*sin(g*k2) - 
	(90*pow(ap,2)*pow(k2,2)*sin(g*k2))/pow(g,2) + 12*ap*bp*pow(k2,3)*sin(g*k2) + 
	7*pow(ap,2)*pow(k2,4)*sin(g*k2) + 2*(pow(cp,2) + 2*bp*dp)*gsl_sf_Si(g*k1) - 
	2*(pow(cp,2) + 2*bp*dp)*gsl_sf_Si(g*k2)))/pow(g,3);
}
*/
// LP Energy Calculation recipes ... Sep 07 (Tue) 2022

