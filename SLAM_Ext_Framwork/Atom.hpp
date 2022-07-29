#ifndef __ATOM
#define __ATOM

#include <Eigen/Core>
#include <iostream>
#include <string>

static int cnt = 0;

class Atom
{

friend class Cell;				// Allowing access previlege to 'Cell' class, i.e., Cell can access privates of Atom
friend class Manager;

private:
	double charge;				// Atom charge

	std::string species;			// Atom name
	std::string type;			// Atom description - core, shel, lone (lonepair)

	Eigen::Vector3d frac;			// fractional coord
	Eigen::Vector3d cart;			// fractional -> Cartesian

	Eigen::Vector3d cart_gd;		// g.d. -> geometric derivative
	Eigen::Vector3d cart_gd_int;		// g.d. -> internal ; same as what GULP gives

public:

	Atom( const double frac_x, const double frac_y, const double frac_z, std::string species, std::string type, const Eigen::Vector3d* lattice_vector )
	{
		frac(0) = frac_x;
		frac(1) = frac_y;
		frac(2) = frac_z;

		this->species = species;
		this->type    = type;

		cart = frac(0)*lattice_vector[0] + frac(1)*lattice_vector[1] + frac(2)*lattice_vector[2];
		// cartesian r = x_frac * a + y_frac * b + z_frac * c;
	}

	virtual void SetFeature( const double core_charge, const double shel_charge=0, const double spring_k2=0, const double spring_k4=0)		// VirtualOverriding - must be in the same format
	{
		this->charge = core_charge;
	}

	virtual void ShowFrac() const
	{
		printf("%4.3s%8.4s%12.6lf%12.6lf%12.6lf%12.6lf\n",species.c_str(),"core",frac(0),frac(1),frac(2),charge);
	}

	virtual void ShowCart() const
	{
		printf("%4.3s%8.4s%12.6lf%12.6lf%12.6lf%12.6lf\n",species.c_str(),"core",cart(0),cart(1),cart(2),charge);
	}

	virtual void InitialiseDerivative()
	{
		this->cart_gd.setZero();
		this->cart_gd_int.setZero();
	}

	virtual void UpdateDerivativeInternal( const Eigen::Matrix3d& lattice_matrix )
	{
		this->cart_gd_int = lattice_matrix * this->cart_gd;
	}
	
	virtual ~Atom()								// Explicit Destructor Check
	{	
		cnt++;
		std::cout << "Atom Destructor ~" << cnt << std::endl;
	}

	std::string GetSpecies() const { return species; }
	std::string GetType() const { return type; }
};

class Shell : public Atom
{

friend class Cell;				// Allowing access previlege to 'Cell' class, i.e., Cell can access privates of Atom
friend class Manager;

private:

	double shel_charge;			// Shell charge
	double spring_k2;
	double spring_k4;

	Eigen::Vector3d shel_frac;		// Below all same but for 'Shell'
	Eigen::Vector3d shel_cart;
	
	Eigen::Vector3d shel_cart_gd;
	Eigen::Vector3d shel_cart_gd_int;

public:

	Shell( const double frac_x, const double frac_y, const double frac_z, std::string species, std::string type, const Eigen::Vector3d* lattice_vector,
	       const double shel_frac_x, const double shel_frac_y, const double shel_frac_z )
	: Atom(frac_x,frac_y,frac_z,species,type,lattice_vector)
	{
		shel_frac(0) = shel_frac_x;
		shel_frac(1) = shel_frac_y;
		shel_frac(2) = shel_frac_z;
		shel_cart = shel_frac(0)*lattice_vector[0] + shel_frac(1)*lattice_vector[1] + shel_frac(2)*lattice_vector[2];
	}

	virtual void SetFeature( const double core_charge, const double shel_charge=0, const double spring_k2=0, const double spring_k4=0 ) override		// VirtualOverriding - must be in the same format
	{																		// Or make pure virutal funcions first in the Base class
		Atom::SetFeature(core_charge);
		this->shel_charge = shel_charge;
		this->spring_k2   = spring_k2;
		this->spring_k4   = spring_k4;
	}

	virtual void ShowFrac() const override
	{	
		Atom::ShowFrac();
		printf("%4.3s%8.4s%12.6lf%12.6lf%12.6lf%12.6lf%12.4lf/%.4lf\n",this->GetSpecies().c_str(),this->GetType().c_str(),this->shel_frac(0),this->shel_frac(1),this->shel_frac(2),
			this->shel_charge,this->spring_k2,this->spring_k4);
	}

	virtual void ShowCart() const override
	{
		Atom::ShowCart();
		printf("%4.3s%8.4s%12.6lf%12.6lf%12.6lf%12.6lf%12.4lf/%.4lf\n",this->GetSpecies().c_str(),this->GetType().c_str(),this->shel_cart(0),this->shel_cart(1),this->shel_cart(2),
			this->shel_charge,this->spring_k2,this->spring_k4);
	}

	virtual void InitialiseDerivative() override
	{	
		Atom::InitialiseDerivative();
		this->shel_cart_gd.setZero();
		this->shel_cart_gd_int.setZero();
	}

	virtual void UpdateDerivativeInternal( const Eigen::Matrix3d& lattice_matrix )
	{
		Atom::UpdateDerivativeInternal( lattice_matrix );
		this->shel_cart_gd_int = lattice_matrix * this->shel_cart_gd;
	}

	virtual ~Shell()
	{
		cnt++;
		std::cout << "Shell Destructor ~" << cnt << std::endl;
	}
};

class LonePair : public Atom
{



};

#endif
