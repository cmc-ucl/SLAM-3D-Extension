#include "Manager.hpp"

#include <iostream>
#include <cmath>
#include <Eigen/Dense>

void Manager::InitialiseEnergy( Cell& C )
{
	C.energy_real_sum_cnt = 0;
	C.energy_reci_sum_cnt = 0;
	C.mono_real_energy = C.mono_reci_energy = C.mono_reci_self_energy = C.mono_total_energy = 0.;
}

void Manager::InitialiseDerivative( Cell& C )
{	
	C.derivative_real_sum_cnt = 0;
	C.derivative_reci_sum_cnt = 0;
	for(int i=0;i<C.NumberOfAtoms;i++) {	C.AtomList[i]->InitialiseDerivative();	}		// Initialise Derivative Field
	C.lattice_sd.setZero();										// Initialise Strain Drivative Field
}

////	Optimise Periodic Summation Workload

void Manager::InitialisePeriodicSysParameter( Cell& C )		// Prepare Parameters - Periodic Summation
{
	int NumberOfObject = 0;

	for(int i=0;i<C.NumberOfAtoms;i++)
	{	NumberOfObject++;
		if( C.AtomList[i]->type == "shel" ) { NumberOfObject++; }
	}

	C.sigma = pow(C.weight*NumberOfObject*M_PI*M_PI*M_PI/C.volume/C.volume,-1./6.);
	C.rcut  = std::sqrt(-log(C.accuracy)*C.sigma*C.sigma);
	C.gcut  = 2./C.sigma*std::sqrt(-log(C.accuracy));


	C.h_max  = static_cast<int>(C.rcut / C.real_vector[0].norm());
	C.k_max  = static_cast<int>(C.rcut / C.real_vector[1].norm());
	C.l_max  = static_cast<int>(C.rcut / C.real_vector[2].norm());
	C.ih_max = static_cast<int>(C.gcut / C.reci_vector[0].norm());
	C.ik_max = static_cast<int>(C.gcut / C.reci_vector[1].norm());
	C.il_max = static_cast<int>(C.gcut / C.reci_vector[2].norm());
}

								// reserved for : Cell.lp_transformation_matrix, AtomList[i/j]->cart, static_cast<Shell*>(AtomList[i/j])->shel_cart,
const Eigen::Matrix4d& Manager::LonePairGetTransformationMatrix( Eigen::Matrix4d& transform_matrix /*in/out*/, const Eigen::Vector3d cart_i, const Eigen::Vector3d cart_j )
{
	//Eigen::Vector3d Rij = cart_i - cart_j;
	Eigen::Vector3d Rij = cart_j - cart_i;			// Using This Convention ... Must follow that 'cart_i' has to be a core position of LonePair of interest
								//					      'cart_j' has to be a classic core / or the other LP
	const double Rxy = sqrt(Rij(0)*Rij(0) + Rij(1)*Rij(1));
	const double R   = Rij.norm();
	double n1, n2, tmp;// dummy variables for workspace

	transform_matrix.setZero();// init Return

	transform_matrix(0,0) = 1.;

	if( (Rij(0) == 0) && (Rij(1) == 0) && (Rij(2) > 0) )			// if vector 'Rij' is on z-axis
	{	transform_matrix(1,1) = 1.;
		transform_matrix(2,2) = 1.;
		transform_matrix(3,3) = 1.;				// set lower-right block 3x3 matrix as I
	}
	else if( (Rij(0) == 0) && (Rij(1) == 0) && (Rij(2) < 0) )		// if vector 'Rij' is on negative z-axis
	{	transform_matrix(1,1) = 1.;
		transform_matrix(2,2) = 1.;
		transform_matrix(3,3) =-1.;				// set the matrix has xy-plane reflection
	}
	else
	{	transform_matrix(3,1) = Rij(0)/R;
		transform_matrix(3,2) = Rij(1)/R;
		transform_matrix(3,3) = Rij(1)/R;			// set local z-axis (k') in the transformed symmetry

		transform_matrix(2,1) = Rij(2)*Rij(0)/Rxy;	
		transform_matrix(2,2) = Rij(2)*Rij(1)/Rxy;	
		transform_matrix(2,3) = -R*sqrt(1.-Rij(2)*Rij(2)/R/R);
	
		n1 = 1./sqrt(transform_matrix(2,1)*transform_matrix(2,1) + transform_matrix(2,2)*transform_matrix(2,2) + transform_matrix(2,3)*transform_matrix(2,3));

		for(int k=0;k<3;k++)
		{	tmp = transform_matrix(2,k+1);
			tmp = tmp*n1;
			transform_matrix(2,k+1) = tmp;			// set local y-axis (j') in the transformed symmetry
		}

		transform_matrix(1,1) = n1/R * ( Rij(2)*Rij(2)*Rij(1)/Rxy + R*Rij(1)*sqrt(1.-Rij(2)*Rij(2)/R/R) );
		transform_matrix(1,2) = n1/R * (-Rij(2)*Rij(2)*Rij(0)/Rxy - R*Rij(0)*sqrt(1.-Rij(2)*Rij(2)/R/R) );
		transform_matrix(1,3) = 0.;

		n2 = 1./sqrt(transform_matrix(1,1)*transform_matrix(1,1) + transform_matrix(1,2)*transform_matrix(1,2) + transform_matrix(1,3)*transform_matrix(1,3));

		for(int k=0;k<3;k++)
		{	tmp = transform_matrix(1,k+1);
			tmp = tmp/n2;
			transform_matrix(1,k+1) = tmp;
		}
	}

	return transform_matrix;
}


////	////	////	////	////	////	////	////	////	////	////	////	////

////	Coulomb Interaction ( Periodic Summation )

void Manager::CoulombMonoMonoReal( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector )
{
	double Qi,Qj;
	Eigen::Vector3d Rij;
        // TransVector = h*a + k*b + l*c
        // Rij         = Ai.r - Aj.r - TransVector;

        if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {       
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;

		C.mono_real_energy += 0.5*(Qi*Qj)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;
        }       

	if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "shel" ) 
        {
		// 1. Handling Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		C.mono_real_energy += 0.5*(Qi*Qj)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;
		// 2. Handling Core - Shel
		Qi = C.AtomList[i]->charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart - TransVector;
		C.mono_real_energy += 0.5*(Qi*Qj)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;
        }       
        
	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "core" ) 
        {
		// 1. Handling Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		C.mono_real_energy += 0.5*(Qi*Qj)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;
		// 2. Handling Shel - Core
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = C.AtomList[j]->charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - C.AtomList[j]->cart - TransVector;
		C.mono_real_energy += 0.5*(Qi*Qj)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;
        }       
        
	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "shel" ) 
        {
		// 1. Handling Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		C.mono_real_energy += 0.5*(Qi*Qj)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;
		// 2. Handling Core - Shel
		Qi = C.AtomList[i]->charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart - TransVector;
		C.mono_real_energy += 0.5*(Qi*Qj)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;
		// 3. Handling Shel - Core
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = C.AtomList[j]->charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - C.AtomList[j]->cart - TransVector;
		C.mono_real_energy += 0.5*(Qi*Qj)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;
		// 4. Handling Shel - Shel
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - static_cast<Shell*>(C.AtomList[j])->shel_cart - TransVector;
		C.mono_real_energy += 0.5*(Qi*Qj)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;
        }       
}       

void Manager::CoulombMonoMonoSelf( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector )
{
	double Qic,Qjc,Qis,Qjs;
	Eigen::Vector3d Rij;
        // TransVector = h*a + k*b + l*c
        // Rij         = Ai.r - Aj.r;
	
        if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {       
		Qic = C.AtomList[i]->charge;
		Qjc = C.AtomList[j]->charge;
		C.mono_reci_self_energy += -0.5*(Qic*Qjc)*2./C.sigma/sqrt(M_PI) * C.TO_EV;
        }

	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "shel" ) 
        {
		// 1. Handling Core - Core
		Qic = C.AtomList[i]->charge;
		Qis = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qjc = C.AtomList[j]->charge;
		Qjs = static_cast<Shell*>(C.AtomList[j])->shel_charge;

		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - static_cast<Shell*>(C.AtomList[j])->cart;

		if ( Rij.norm() != 0. )	// if shell / core is seperated !
		{
			C.mono_reci_self_energy += -0.5*(Qic*Qjc + Qis*Qis)*2./C.sigma/sqrt(M_PI) * C.TO_EV;
			C.mono_reci_self_energy += -(Qic*Qjs)/Rij.norm()*erf(Rij.norm()/C.sigma)* C.TO_EV;
		}
		else
		{
			C.mono_reci_self_energy += -0.5*(Qic*Qjc + Qic*Qis + Qis*Qic + Qis*Qis)*2./C.sigma/sqrt(M_PI) * C.TO_EV;
		}
			
        }       
}       

void Manager::CoulombMonoMonoReci( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector )
{
	double Qi,Qj;
	double g_norm = TransVector.norm();
	double g_sqr  = g_norm*g_norm;
	Eigen::Vector3d Rij;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r;
        
        if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {	
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;
		C.mono_reci_energy += (2.*M_PI/C.volume)*(Qi*Qj)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;
        }       
        
	if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "shel" ) 
        {
		// 1. Handling Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;
		C.mono_reci_energy += (2.*M_PI/C.volume)*(Qi*Qj)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;
		// 2. Handling Core - Shel
		Qi = C.AtomList[i]->charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart;
		C.mono_reci_energy += (2.*M_PI/C.volume)*(Qi*Qj)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;
        }       
        
	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "core" ) 
        {
		// 1. Handling Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;
		C.mono_reci_energy += (2.*M_PI/C.volume)*(Qi*Qj)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;
		// 3. Handling Shel - Core
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = C.AtomList[j]->charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - C.AtomList[j]->cart;
		C.mono_reci_energy += (2.*M_PI/C.volume)*(Qi*Qj)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;
        }       
        
	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "shel" ) 
        {
		// 1. Handling Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;
		C.mono_reci_energy += (2.*M_PI/C.volume)*(Qi*Qj)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;
		// 2. Handling Core - Shel
		Qi = C.AtomList[i]->charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart;
		C.mono_reci_energy += (2.*M_PI/C.volume)*(Qi*Qj)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;
		// 3. Handling Shel - Core
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = C.AtomList[j]->charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - C.AtomList[j]->cart;
		C.mono_reci_energy += (2.*M_PI/C.volume)*(Qi*Qj)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;
		// 4. Handling Shel - Shel
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - static_cast<Shell*>(C.AtomList[j])->shel_cart;
		C.mono_reci_energy += (2.*M_PI/C.volume)*(Qi*Qj)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;
        }       
} 

////	////	////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////

////	Geometric (RAW) Derivatives

void Manager::CoulombDerivativeReal( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector )
{
	double Qi,Qj;
	double r_norm,r_sqr;
	Eigen::Vector3d Rij;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r - TransVector;
	double intact;

        if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.AtomList[i]->cart_gd += intact*Rij;
		C.AtomList[j]->cart_gd -= intact*Rij;
        }       
        
	if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "shel" ) 
        {
		// Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;
	
		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.AtomList[i]->cart_gd += intact*Rij;
		C.AtomList[j]->cart_gd -= intact*Rij;
		
		// Core - Shell
		Qi = C.AtomList[i]->charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;
// r_norm, r_sqr were missed!!!!
	
		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.AtomList[i]->cart_gd += intact*Rij;
		static_cast<Shell*>(C.AtomList[j])->shel_cart_gd -= intact*Rij;

        }       
        
	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "core" ) 
        {
		// Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;
	
		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.AtomList[i]->cart_gd += intact*Rij;
		C.AtomList[j]->cart_gd -= intact*Rij;

		// Shel - Core
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = C.AtomList[j]->charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		static_cast<Shell*>(C.AtomList[i])->shel_cart_gd += intact*Rij;
		C.AtomList[j]->cart_gd -= intact*Rij;
        }       
        
	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "shel" ) 
        {
		// Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;
	
		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.AtomList[i]->cart_gd += intact*Rij;
		C.AtomList[j]->cart_gd -= intact*Rij;
		
		// Core - Shell
		Qi = C.AtomList[i]->charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;
// r_norm, r_sqr were missed!!!!

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.AtomList[i]->cart_gd += intact*Rij;
		static_cast<Shell*>(C.AtomList[j])->shel_cart_gd -= intact*Rij;

		// Shell - Core
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = C.AtomList[j]->charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		static_cast<Shell*>(C.AtomList[i])->shel_cart_gd += intact*Rij;
		C.AtomList[j]->cart_gd -= intact*Rij;

		// Shell - Shell
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - static_cast<Shell*>(C.AtomList[j])->shel_cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		static_cast<Shell*>(C.AtomList[i])->shel_cart_gd += intact*Rij;
		static_cast<Shell*>(C.AtomList[j])->shel_cart_gd -= intact*Rij;
        }       
}

void Manager::CoulombDerivativeSelf( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector )
{
	double Qi,Qj;
	double r_norm,r_sqr;
	Eigen::Vector3d Rij;
	
	double intact;

	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "shel" ) 
        {
		// Self - Core/Shell
		Qi = C.AtomList[i]->charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		if ( r_norm != 0. )	// Case of Shell-Core separation is not zero
		{
			intact = 2.0 * (-0.5*Qi*Qj*(2./C.sigma/sqrt(M_PI)*exp(-r_sqr/C.sigma/C.sigma)/r_sqr - erf(r_norm/C.sigma)/r_norm/r_sqr) * C.TO_EV);
			// shell - core , core - shell counting twice

			C.AtomList[i]->cart_gd += intact * Rij;
			static_cast<Shell*>(C.AtomList[j])->shel_cart_gd -= intact * Rij;
		}
	}
}       

void Manager::CoulombDerivativeReci( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector )
{
	double Qi,Qj;
	double g_norm = TransVector.norm();
	double g_sqr  = g_norm*g_norm;
	Eigen::Vector3d Rij;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r;
	double intact;
        
        if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;

		intact = C.TO_EV*((2.*M_PI)/C.volume)*Qi*Qj*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);

		C.AtomList[i]->cart_gd += intact * TransVector;
		C.AtomList[j]->cart_gd -= intact * TransVector;
        }       
        
	if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "shel" ) 
        {	
		// Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;

		intact = C.TO_EV*((2.*M_PI)/C.volume)*Qi*Qj*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);

		C.AtomList[i]->cart_gd += intact * TransVector;
		C.AtomList[j]->cart_gd -= intact * TransVector;

		// Core - Shell
		Qi = C.AtomList[i]->charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart;

		intact = C.TO_EV*((2.*M_PI)/C.volume)*Qi*Qj*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);

		C.AtomList[i]->cart_gd += intact * TransVector;
		static_cast<Shell*>(C.AtomList[j])->shel_cart_gd -= intact * TransVector;
        }       
        
	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "core" ) 
        {
		// Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;

		intact = C.TO_EV*((2.*M_PI)/C.volume)*Qi*Qj*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);

		C.AtomList[i]->cart_gd += intact * TransVector;
		C.AtomList[j]->cart_gd -= intact * TransVector;

		// Shell - Core
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = C.AtomList[j]->charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - C.AtomList[j]->cart;
	
		intact = C.TO_EV*((2.*M_PI)/C.volume)*Qi*Qj*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);

		static_cast<Shell*>(C.AtomList[i])->shel_cart_gd += intact * TransVector;
		C.AtomList[j]->cart_gd -= intact * TransVector;
        }       
        
	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "shel" ) 
        {
		// Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;

		intact = C.TO_EV*((2.*M_PI)/C.volume)*Qi*Qj*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);

		C.AtomList[i]->cart_gd += intact * TransVector;
		C.AtomList[j]->cart_gd -= intact * TransVector;

		// Core - Shell
		Qi = C.AtomList[i]->charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart;

		intact = C.TO_EV*((2.*M_PI)/C.volume)*Qi*Qj*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);

		C.AtomList[i]->cart_gd += intact * TransVector;
		static_cast<Shell*>(C.AtomList[j])->shel_cart_gd -= intact * TransVector;

		// Shell - Core
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = C.AtomList[j]->charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - C.AtomList[j]->cart;
	
		intact = C.TO_EV*((2.*M_PI)/C.volume)*Qi*Qj*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);

		static_cast<Shell*>(C.AtomList[i])->shel_cart_gd += intact * TransVector;
		C.AtomList[j]->cart_gd -= intact * TransVector;

		// Shell - Shell
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - static_cast<Shell*>(C.AtomList[j])->shel_cart;
	
		intact = C.TO_EV*((2.*M_PI)/C.volume)*Qi*Qj*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);

		static_cast<Shell*>(C.AtomList[i])->shel_cart_gd += intact * TransVector;
		static_cast<Shell*>(C.AtomList[j])->shel_cart_gd -= intact * TransVector;
        }       
}       

////	////	////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////

////	Strain Derivative

void Manager::StrainDerivativeReal( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector )
{
	double Qi,Qj;
	double r_norm,r_sqr;
	Eigen::Vector3d Rij;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r - TransVector;
	double intact;
        
        if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.lattice_sd(0,0) += intact * Rij(0) * Rij(0);	C.lattice_sd(0,1) += intact * Rij(0) * Rij(1);	C.lattice_sd(0,2) += intact * Rij(0) * Rij(2);
								C.lattice_sd(1,1) += intact * Rij(1) * Rij(1);	C.lattice_sd(1,2) += intact * Rij(1) * Rij(2);
														C.lattice_sd(2,2) += intact * Rij(2) * Rij(2);
        }       
        
	if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "shel" ) 
        {	
		// Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.lattice_sd(0,0) += intact * Rij(0) * Rij(0);	C.lattice_sd(0,1) += intact * Rij(0) * Rij(1);	C.lattice_sd(0,2) += intact * Rij(0) * Rij(2);
								C.lattice_sd(1,1) += intact * Rij(1) * Rij(1);	C.lattice_sd(1,2) += intact * Rij(1) * Rij(2);
														C.lattice_sd(2,2) += intact * Rij(2) * Rij(2);
		// Core - Shell
		Qi = C.AtomList[i]->charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.lattice_sd(0,0) += intact * Rij(0) * Rij(0);	C.lattice_sd(0,1) += intact * Rij(0) * Rij(1);	C.lattice_sd(0,2) += intact * Rij(0) * Rij(2);
								C.lattice_sd(1,1) += intact * Rij(1) * Rij(1);	C.lattice_sd(1,2) += intact * Rij(1) * Rij(2);
														C.lattice_sd(2,2) += intact * Rij(2) * Rij(2);
        }       
        
	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "core" ) 
        {	
		// Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.lattice_sd(0,0) += intact * Rij(0) * Rij(0);	C.lattice_sd(0,1) += intact * Rij(0) * Rij(1);	C.lattice_sd(0,2) += intact * Rij(0) * Rij(2);
								C.lattice_sd(1,1) += intact * Rij(1) * Rij(1);	C.lattice_sd(1,2) += intact * Rij(1) * Rij(2);
														C.lattice_sd(2,2) += intact * Rij(2) * Rij(2);
		// Shell - Core
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = C.AtomList[j]->charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.lattice_sd(0,0) += intact * Rij(0) * Rij(0);	C.lattice_sd(0,1) += intact * Rij(0) * Rij(1);	C.lattice_sd(0,2) += intact * Rij(0) * Rij(2);
								C.lattice_sd(1,1) += intact * Rij(1) * Rij(1);	C.lattice_sd(1,2) += intact * Rij(1) * Rij(2);
														C.lattice_sd(2,2) += intact * Rij(2) * Rij(2);
        }       
        
	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "shel" ) 
        {
		// Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.lattice_sd(0,0) += intact * Rij(0) * Rij(0);	C.lattice_sd(0,1) += intact * Rij(0) * Rij(1);	C.lattice_sd(0,2) += intact * Rij(0) * Rij(2);
								C.lattice_sd(1,1) += intact * Rij(1) * Rij(1);	C.lattice_sd(1,2) += intact * Rij(1) * Rij(2);
														C.lattice_sd(2,2) += intact * Rij(2) * Rij(2);
		// Core - Shell
		Qi = C.AtomList[i]->charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.lattice_sd(0,0) += intact * Rij(0) * Rij(0);	C.lattice_sd(0,1) += intact * Rij(0) * Rij(1);	C.lattice_sd(0,2) += intact * Rij(0) * Rij(2);
								C.lattice_sd(1,1) += intact * Rij(1) * Rij(1);	C.lattice_sd(1,2) += intact * Rij(1) * Rij(2);
														C.lattice_sd(2,2) += intact * Rij(2) * Rij(2);
		// Shell - Core
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = C.AtomList[j]->charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.lattice_sd(0,0) += intact * Rij(0) * Rij(0);	C.lattice_sd(0,1) += intact * Rij(0) * Rij(1);	C.lattice_sd(0,2) += intact * Rij(0) * Rij(2);
								C.lattice_sd(1,1) += intact * Rij(1) * Rij(1);	C.lattice_sd(1,2) += intact * Rij(1) * Rij(2);
														C.lattice_sd(2,2) += intact * Rij(2) * Rij(2);
		// Shell - Shell
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - static_cast<Shell*>(C.AtomList[j])->shel_cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.lattice_sd(0,0) += intact * Rij(0) * Rij(0);	C.lattice_sd(0,1) += intact * Rij(0) * Rij(1);	C.lattice_sd(0,2) += intact * Rij(0) * Rij(2);
								C.lattice_sd(1,1) += intact * Rij(1) * Rij(1);	C.lattice_sd(1,2) += intact * Rij(1) * Rij(2);
														C.lattice_sd(2,2) += intact * Rij(2) * Rij(2);
        }       
}       

void Manager::StrainDerivativeSelf( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector )
{
	double Qi,Qj;
	double r_norm,r_sqr;
	Eigen::Vector3d Rij;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r - TransVector;
	double intact;

        if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "shel" )    // Handling Core - Core (i.e., charge charge interaction);
	{
		// Self - Core/Shell
		Qi = C.AtomList[i]->charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		if ( r_norm != 0. )	// Case of Shell-Core separation is not zero
		{
			intact = 2.0 * (-0.5*Qi*Qj*(2./C.sigma/sqrt(M_PI)*exp(-r_sqr/C.sigma/C.sigma)/r_sqr - erf(r_norm/C.sigma)/r_norm/r_sqr) * C.TO_EV);
			// leading "2.0 *" shell - core , core - shell counting twice

			C.lattice_sd(0,0) += intact * Rij(0) * Rij(0);	C.lattice_sd(0,1) += intact * Rij(0) * Rij(1);	C.lattice_sd(0,2) += intact * Rij(0) * Rij(2);
									C.lattice_sd(1,1) += intact * Rij(1) * Rij(1);	C.lattice_sd(1,2) += intact * Rij(1) * Rij(2);
															C.lattice_sd(2,2) += intact * Rij(2) * Rij(2);
		}

	}
}

void Manager::StrainDerivativeReci( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector )
{
	double Qi,Qj;
	double g_norm = TransVector.norm();
	double g_sqr  = g_norm*g_norm;
	Eigen::Vector3d Rij;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r;
	double intact[4];

        if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;

		intact[0] = C.TO_EV*((2.*M_PI)/C.volume)*(Qi*Qj);
		intact[1] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);
		intact[2] = (-2.*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr/g_sqr*cos(TransVector.adjoint()*Rij)-0.5*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*C.sigma*C.sigma*cos(TransVector.adjoint()*Rij));
		intact[3] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*sin(TransVector.adjoint()*Rij);
		intact[4] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*cos(TransVector.adjoint()*Rij);

		// Strain derivative (1) - w.r.t. r_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*intact[1] * TransVector(0) * Rij(0);	C.lattice_sd(0,1) += intact[0]*intact[1] * TransVector(0) * Rij(1);	C.lattice_sd(0,2) += intact[0]*intact[1] * TransVector(0) * Rij(2);
											C.lattice_sd(1,1) += intact[0]*intact[1] * TransVector(1) * Rij(1);	C.lattice_sd(1,2) += intact[0]*intact[1] * TransVector(1) * Rij(2);
																				C.lattice_sd(2,2) += intact[0]*intact[1] * TransVector(2) * Rij(2);
		// Strain derivative (2) - w.r.t. g_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*(intact[2]*TransVector(0)-intact[3]*Rij(0))*-TransVector(0);	C.lattice_sd(0,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(0);	C.lattice_sd(0,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(0);
														C.lattice_sd(1,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(1);	C.lattice_sd(1,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(1);
																										C.lattice_sd(2,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(2);
		// Strain derivative (3) - w.r.t cell volume in the reciprocal space
		C.lattice_sd(0,0) += -intact[0]*intact[4];	C.lattice_sd(1,1) += -intact[0]*intact[4];	C.lattice_sd(2,2) += -intact[0]*intact[4];
        }       
        
	if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "shel" ) 
        {
		// Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;

		intact[0] = C.TO_EV*((2.*M_PI)/C.volume)*(Qi*Qj);
		intact[1] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);
		intact[2] = (-2.*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr/g_sqr*cos(TransVector.adjoint()*Rij)-0.5*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*C.sigma*C.sigma*cos(TransVector.adjoint()*Rij));
		intact[3] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*sin(TransVector.adjoint()*Rij);
		intact[4] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*cos(TransVector.adjoint()*Rij);

		// Strain derivative (1) - w.r.t. r_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*intact[1] * TransVector(0) * Rij(0);	C.lattice_sd(0,1) += intact[0]*intact[1] * TransVector(0) * Rij(1);	C.lattice_sd(0,2) += intact[0]*intact[1] * TransVector(0) * Rij(2);
											C.lattice_sd(1,1) += intact[0]*intact[1] * TransVector(1) * Rij(1);	C.lattice_sd(1,2) += intact[0]*intact[1] * TransVector(1) * Rij(2);
																				C.lattice_sd(2,2) += intact[0]*intact[1] * TransVector(2) * Rij(2);
		// Strain derivative (2) - w.r.t. g_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*(intact[2]*TransVector(0)-intact[3]*Rij(0))*-TransVector(0);	C.lattice_sd(0,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(0);	C.lattice_sd(0,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(0);
														C.lattice_sd(1,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(1);	C.lattice_sd(1,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(1);
																										C.lattice_sd(2,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(2);
		// Strain derivative (3) - w.r.t cell volume in the reciprocal space
		C.lattice_sd(0,0) += -intact[0]*intact[4];	C.lattice_sd(1,1) += -intact[0]*intact[4];	C.lattice_sd(2,2) += -intact[0]*intact[4];

		// Core - Shell
		Qi = C.AtomList[i]->charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart;

		intact[0] = C.TO_EV*((2.*M_PI)/C.volume)*(Qi*Qj);
		intact[1] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);
		intact[2] = (-2.*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr/g_sqr*cos(TransVector.adjoint()*Rij)-0.5*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*C.sigma*C.sigma*cos(TransVector.adjoint()*Rij));
		intact[3] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*sin(TransVector.adjoint()*Rij);
		intact[4] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*cos(TransVector.adjoint()*Rij);

		// Strain derivative (1) - w.r.t. r_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*intact[1] * TransVector(0) * Rij(0);	C.lattice_sd(0,1) += intact[0]*intact[1] * TransVector(0) * Rij(1);	C.lattice_sd(0,2) += intact[0]*intact[1] * TransVector(0) * Rij(2);
											C.lattice_sd(1,1) += intact[0]*intact[1] * TransVector(1) * Rij(1);	C.lattice_sd(1,2) += intact[0]*intact[1] * TransVector(1) * Rij(2);
																				C.lattice_sd(2,2) += intact[0]*intact[1] * TransVector(2) * Rij(2);
		// Strain derivative (2) - w.r.t. g_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*(intact[2]*TransVector(0)-intact[3]*Rij(0))*-TransVector(0);	C.lattice_sd(0,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(0);	C.lattice_sd(0,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(0);
														C.lattice_sd(1,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(1);	C.lattice_sd(1,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(1);
																										C.lattice_sd(2,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(2);
		// Strain derivative (3) - w.r.t cell volume in the reciprocal space
		C.lattice_sd(0,0) += -intact[0]*intact[4];	C.lattice_sd(1,1) += -intact[0]*intact[4];	C.lattice_sd(2,2) += -intact[0]*intact[4];
        }       
        
	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "core" ) 
        {
		// Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;

		intact[0] = C.TO_EV*((2.*M_PI)/C.volume)*(Qi*Qj);
		intact[1] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);
		intact[2] = (-2.*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr/g_sqr*cos(TransVector.adjoint()*Rij)-0.5*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*C.sigma*C.sigma*cos(TransVector.adjoint()*Rij));
		intact[3] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*sin(TransVector.adjoint()*Rij);
		intact[4] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*cos(TransVector.adjoint()*Rij);

		// Strain derivative (1) - w.r.t. r_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*intact[1] * TransVector(0) * Rij(0);	C.lattice_sd(0,1) += intact[0]*intact[1] * TransVector(0) * Rij(1);	C.lattice_sd(0,2) += intact[0]*intact[1] * TransVector(0) * Rij(2);
											C.lattice_sd(1,1) += intact[0]*intact[1] * TransVector(1) * Rij(1);	C.lattice_sd(1,2) += intact[0]*intact[1] * TransVector(1) * Rij(2);
																				C.lattice_sd(2,2) += intact[0]*intact[1] * TransVector(2) * Rij(2);
		// Strain derivative (2) - w.r.t. g_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*(intact[2]*TransVector(0)-intact[3]*Rij(0))*-TransVector(0);	C.lattice_sd(0,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(0);	C.lattice_sd(0,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(0);
														C.lattice_sd(1,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(1);	C.lattice_sd(1,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(1);
																										C.lattice_sd(2,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(2);
		// Strain derivative (3) - w.r.t cell volume in the reciprocal space
		C.lattice_sd(0,0) += -intact[0]*intact[4];	C.lattice_sd(1,1) += -intact[0]*intact[4];	C.lattice_sd(2,2) += -intact[0]*intact[4];

		// Shell - Core
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = C.AtomList[j]->charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - C.AtomList[j]->cart;

		intact[0] = C.TO_EV*((2.*M_PI)/C.volume)*(Qi*Qj);
		intact[1] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);
		intact[2] = (-2.*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr/g_sqr*cos(TransVector.adjoint()*Rij)-0.5*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*C.sigma*C.sigma*cos(TransVector.adjoint()*Rij));
		intact[3] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*sin(TransVector.adjoint()*Rij);
		intact[4] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*cos(TransVector.adjoint()*Rij);

		// Strain derivative (1) - w.r.t. r_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*intact[1] * TransVector(0) * Rij(0);	C.lattice_sd(0,1) += intact[0]*intact[1] * TransVector(0) * Rij(1);	C.lattice_sd(0,2) += intact[0]*intact[1] * TransVector(0) * Rij(2);
											C.lattice_sd(1,1) += intact[0]*intact[1] * TransVector(1) * Rij(1);	C.lattice_sd(1,2) += intact[0]*intact[1] * TransVector(1) * Rij(2);
																				C.lattice_sd(2,2) += intact[0]*intact[1] * TransVector(2) * Rij(2);
		// Strain derivative (2) - w.r.t. g_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*(intact[2]*TransVector(0)-intact[3]*Rij(0))*-TransVector(0);	C.lattice_sd(0,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(0);	C.lattice_sd(0,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(0);
														C.lattice_sd(1,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(1);	C.lattice_sd(1,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(1);
																										C.lattice_sd(2,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(2);
		// Strain derivative (3) - w.r.t cell volume in the reciprocal space
		C.lattice_sd(0,0) += -intact[0]*intact[4];	C.lattice_sd(1,1) += -intact[0]*intact[4];	C.lattice_sd(2,2) += -intact[0]*intact[4];
        }       
        
	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "shel" ) 
        {
		// Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;

		intact[0] = C.TO_EV*((2.*M_PI)/C.volume)*(Qi*Qj);
		intact[1] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);
		intact[2] = (-2.*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr/g_sqr*cos(TransVector.adjoint()*Rij)-0.5*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*C.sigma*C.sigma*cos(TransVector.adjoint()*Rij));
		intact[3] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*sin(TransVector.adjoint()*Rij);
		intact[4] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*cos(TransVector.adjoint()*Rij);

		// Strain derivative (1) - w.r.t. r_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*intact[1] * TransVector(0) * Rij(0);	C.lattice_sd(0,1) += intact[0]*intact[1] * TransVector(0) * Rij(1);	C.lattice_sd(0,2) += intact[0]*intact[1] * TransVector(0) * Rij(2);
											C.lattice_sd(1,1) += intact[0]*intact[1] * TransVector(1) * Rij(1);	C.lattice_sd(1,2) += intact[0]*intact[1] * TransVector(1) * Rij(2);
																				C.lattice_sd(2,2) += intact[0]*intact[1] * TransVector(2) * Rij(2);
		// Strain derivative (2) - w.r.t. g_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*(intact[2]*TransVector(0)-intact[3]*Rij(0))*-TransVector(0);	C.lattice_sd(0,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(0);	C.lattice_sd(0,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(0);
														C.lattice_sd(1,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(1);	C.lattice_sd(1,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(1);
																										C.lattice_sd(2,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(2);
		// Strain derivative (3) - w.r.t cell volume in the reciprocal space
		C.lattice_sd(0,0) += -intact[0]*intact[4];	C.lattice_sd(1,1) += -intact[0]*intact[4];	C.lattice_sd(2,2) += -intact[0]*intact[4];

		// Core - Shell
		Qi = C.AtomList[i]->charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart;

		intact[0] = C.TO_EV*((2.*M_PI)/C.volume)*(Qi*Qj);
		intact[1] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);
		intact[2] = (-2.*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr/g_sqr*cos(TransVector.adjoint()*Rij)-0.5*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*C.sigma*C.sigma*cos(TransVector.adjoint()*Rij));
		intact[3] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*sin(TransVector.adjoint()*Rij);
		intact[4] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*cos(TransVector.adjoint()*Rij);

		// Strain derivative (1) - w.r.t. r_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*intact[1] * TransVector(0) * Rij(0);	C.lattice_sd(0,1) += intact[0]*intact[1] * TransVector(0) * Rij(1);	C.lattice_sd(0,2) += intact[0]*intact[1] * TransVector(0) * Rij(2);
											C.lattice_sd(1,1) += intact[0]*intact[1] * TransVector(1) * Rij(1);	C.lattice_sd(1,2) += intact[0]*intact[1] * TransVector(1) * Rij(2);
																				C.lattice_sd(2,2) += intact[0]*intact[1] * TransVector(2) * Rij(2);
		// Strain derivative (2) - w.r.t. g_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*(intact[2]*TransVector(0)-intact[3]*Rij(0))*-TransVector(0);	C.lattice_sd(0,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(0);	C.lattice_sd(0,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(0);
														C.lattice_sd(1,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(1);	C.lattice_sd(1,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(1);
																										C.lattice_sd(2,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(2);
		// Strain derivative (3) - w.r.t cell volume in the reciprocal space
		C.lattice_sd(0,0) += -intact[0]*intact[4];	C.lattice_sd(1,1) += -intact[0]*intact[4];	C.lattice_sd(2,2) += -intact[0]*intact[4];

		// Shell - Core
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = C.AtomList[j]->charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - C.AtomList[j]->cart;

		intact[0] = C.TO_EV*((2.*M_PI)/C.volume)*(Qi*Qj);
		intact[1] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);
		intact[2] = (-2.*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr/g_sqr*cos(TransVector.adjoint()*Rij)-0.5*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*C.sigma*C.sigma*cos(TransVector.adjoint()*Rij));
		intact[3] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*sin(TransVector.adjoint()*Rij);
		intact[4] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*cos(TransVector.adjoint()*Rij);

		// Strain derivative (1) - w.r.t. r_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*intact[1] * TransVector(0) * Rij(0);	C.lattice_sd(0,1) += intact[0]*intact[1] * TransVector(0) * Rij(1);	C.lattice_sd(0,2) += intact[0]*intact[1] * TransVector(0) * Rij(2);
											C.lattice_sd(1,1) += intact[0]*intact[1] * TransVector(1) * Rij(1);	C.lattice_sd(1,2) += intact[0]*intact[1] * TransVector(1) * Rij(2);
																				C.lattice_sd(2,2) += intact[0]*intact[1] * TransVector(2) * Rij(2);
		// Strain derivative (2) - w.r.t. g_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*(intact[2]*TransVector(0)-intact[3]*Rij(0))*-TransVector(0);	C.lattice_sd(0,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(0);	C.lattice_sd(0,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(0);
														C.lattice_sd(1,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(1);	C.lattice_sd(1,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(1);
																										C.lattice_sd(2,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(2);
		// Strain derivative (3) - w.r.t cell volume in the reciprocal space
		C.lattice_sd(0,0) += -intact[0]*intact[4];	C.lattice_sd(1,1) += -intact[0]*intact[4];	C.lattice_sd(2,2) += -intact[0]*intact[4];
		
		// Shell - Shell
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - static_cast<Shell*>(C.AtomList[j])->shel_cart;

		intact[0] = C.TO_EV*((2.*M_PI)/C.volume)*(Qi*Qj);
		intact[1] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);
		intact[2] = (-2.*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr/g_sqr*cos(TransVector.adjoint()*Rij)-0.5*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*C.sigma*C.sigma*cos(TransVector.adjoint()*Rij));
		intact[3] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*sin(TransVector.adjoint()*Rij);
		intact[4] = exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr*cos(TransVector.adjoint()*Rij);

		// Strain derivative (1) - w.r.t. r_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*intact[1] * TransVector(0) * Rij(0);	C.lattice_sd(0,1) += intact[0]*intact[1] * TransVector(0) * Rij(1);	C.lattice_sd(0,2) += intact[0]*intact[1] * TransVector(0) * Rij(2);
											C.lattice_sd(1,1) += intact[0]*intact[1] * TransVector(1) * Rij(1);	C.lattice_sd(1,2) += intact[0]*intact[1] * TransVector(1) * Rij(2);
																				C.lattice_sd(2,2) += intact[0]*intact[1] * TransVector(2) * Rij(2);
		// Strain derivative (2) - w.r.t. g_vector in the reciprocal space
		C.lattice_sd(0,0) += intact[0]*(intact[2]*TransVector(0)-intact[3]*Rij(0))*-TransVector(0);	C.lattice_sd(0,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(0);	C.lattice_sd(0,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(0);
														C.lattice_sd(1,1) += intact[0]*(intact[2]*TransVector(1)-intact[3]*Rij(1))*-TransVector(1);	C.lattice_sd(1,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(1);
																										C.lattice_sd(2,2) += intact[0]*(intact[2]*TransVector(2)-intact[3]*Rij(2))*-TransVector(2);
		// Strain derivative (3) - w.r.t cell volume in the reciprocal space
		C.lattice_sd(0,0) += -intact[0]*intact[4];	C.lattice_sd(1,1) += -intact[0]*intact[4];	C.lattice_sd(2,2) += -intact[0]*intact[4];
        }       
}       
