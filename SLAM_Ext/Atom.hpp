#ifndef __ATOM
#define __ATOM

#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include <sstream>
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
	// READ ONLY METHODS
	const double& GetCoreCharge() const { return this->charge; }
	const std::string& GetSpecies() const { return species; }
	const std::string& GetType() const { return type; }

	const Eigen::Vector3d& GetCoreCart() const { return this->cart; }
	const Eigen::Vector3d& GetCoreFrac() const { return this->frac; }
	const Eigen::Vector3d& GetCoreCartDerivative() const { return this->cart_gd; }
	const Eigen::Vector3d& GetCoreCoartDerivativeInternal() const { return this->cart_gd_int; }

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

	virtual void SetFeature( const double core_charge, const double extra_charge=0, const double val1=0, const double val2=0)		// VirtualOverriding - must be in the same format
	{
		this->charge = core_charge;
	}

	virtual void ShowFrac() const
	{
		printf("%4.3s%8.4s%12.6lf%12.6lf%12.6lf%8.3lf\n",species.c_str(),"core",frac(0),frac(1),frac(2),charge);
	}

	virtual void ShowCart() const
	{
		printf("%4.3s%8.4s%12.6lf%12.6lf%12.6lf%8.3lf\n",species.c_str(),"core",cart(0),cart(1),cart(2),charge);
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

	// Virtuals for later uses
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
	// READ ONLY METHODS
	const double& GetShellCharge() const { return this->shel_charge; }
	const double& GetShellSpringK2() const { return this->spring_k2; }
	const double& GetShellSpringK4() const { return this->spring_k4; }

	const Eigen::Vector3d& GetShellCart() const { return this->shel_cart; }
	const Eigen::Vector3d& GetShellFrac() const { return this->shel_frac; }
	const Eigen::Vector3d& GetShellCartDerivative() const { return this->shel_cart_gd; }
	const Eigen::Vector3d& GetShellCoartDerivativeInternal() const { return this->shel_cart_gd_int; }

	Shell( const double frac_x, const double frac_y, const double frac_z, std::string species, std::string type, const Eigen::Vector3d* lattice_vector,
	       const double shel_frac_x, const double shel_frac_y, const double shel_frac_z )
	: Atom(frac_x,frac_y,frac_z,species,type,lattice_vector)
	{
		shel_frac(0) = shel_frac_x;
		shel_frac(1) = shel_frac_y;
		shel_frac(2) = shel_frac_z;
		shel_cart = shel_frac(0)*lattice_vector[0] + shel_frac(1)*lattice_vector[1] + shel_frac(2)*lattice_vector[2];
	}

	virtual void SetFeature( const double core_charge, const double extra_charge=0, const double val1=0, const double val2=0 ) override		// VirtualOverriding - must be in the same format
	{																		// Or make pure virutal funcions first in the Base class
		Atom::SetFeature(core_charge);
		this->shel_charge = extra_charge;
		this->spring_k2   = val1;
		this->spring_k4   = val2;
	}

	virtual void ShowFrac() const override
	{	
		Atom::ShowFrac();
		printf("%4.3s%8.4s%12.6lf%12.6lf%12.6lf%8.3lf",this->GetSpecies().c_str(),this->GetType().c_str(),this->shel_frac(0),this->shel_frac(1),this->shel_frac(2),this->shel_charge);
		printf("%12.3lf%10.3lf%12.10s%14.13s\n",this->spring_k2,this->spring_k4,"k2(eV/A^2)","/ k4(eV/A^4)");
	}

	virtual void ShowCart() const override
	{
		Atom::ShowCart();
		printf("%4.3s%8.4s%12.6lf%12.6lf%12.6lf%8.3lf",this->GetSpecies().c_str(),this->GetType().c_str(),this->shel_cart(0),this->shel_cart(1),this->shel_cart(2),this->shel_charge);
		printf("%12.3lf%10.3lf%12.10s%14.13s\n",this->spring_k2,this->spring_k4,"k2(eV/A^2)","/ k4(eV/A^4)");
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


#define RADIAL_Pb_KNOT "lp_radial/Pb_radial_knot.dat"
#define RADIAL_Pb_S    "lp_radial/Pb_radial_s.dat"
#define RADIAL_Pb_P    "lp_radial/Pb_radial_p.dat"

#define RADIAL_Sn_KNOT "lp_radial/Sn_radial_knot.dat"
#define RADIAL_Sn_S    "lp_radial/Sn_radial_s.dat"
#define RADIAL_Sn_P    "lp_radial/Sn_radial_p.dat"

#define RADIAL_Bi_KNOT "lp_radial/Bi_radial_knot.dat"
#define RADIAL_Bi_S    "lp_radial/Bi_radial_s.dat"
#define RADIAL_Bi_P    "lp_radial/Bi_radial_p.dat"

class LonePair : public Atom
{

friend class Cell;				// Allowing access previlege to 'Cell' class, i.e., Cell can access privates of Atom
friend class Manager;
//friend class LonePairMatrix_H;

private:

	double lp_charge;	// Lone pair charge
	double lp_lambda;	// Lone pair s <-> p energy sepration

	int lp_gs_index;		// Ground State Index		// This variable is going to be set by manager class member function
	Eigen::EigenSolver<Eigen::Matrix4d> lp_eigensolver;	// EigenSolver itself contains evals/evecs ... accessed by .eigenvalues() / .eigenvectors() // elements of vectors and matrices with type of std::complex<double?>
	/*
		! Solve EigenSystem       : this->lp_eigensolver.compute( $TARGE_MATRIX(SYSTEM) , compute_eigenvectors = true );

		! Access Computed Results : this->lp_eigensolver.eigenvalues()[i].real() ...    real part of the 'i'th eigenvalue

		!			    this->lp_eigensolver.eigenvectors()(i,j).real() ... real part of the (i,j)th eigenvector

// Below should be implemented as the class internal method
cout << "Get Min index" << endl;
std::vector<double> evals;
int min_idx;
for(int i=0;i<3;i++)
{       evals.push_back(es.eigenvalues()(i).real());
}
min_idx = std::min_element(evals.begin(),evals.end()) - evals.begin();
cout << min_idx << endl;
	*/

	Eigen::Matrix4d lp_h_matrix;			// Lone pair Hamiltonian Matrix;
	Eigen::Matrix4d lp_h_matrix_tmp;
	Eigen::Matrix4d lp_h_matrix_derivatives[3];	// 1,2,3 for 'x', 'y', 'z';
	
	Eigen::Vector3d lp_gd;
	Eigen::Vector3d lp_gd_int;			// Gradient on the core by LonePair Density

	// RadialWaveFunctions
	int lp_spline_knot_count;
	std::vector<double> lp_r;
	std::vector<double> lp_r_s_function[4];		// [0-3] : a,b,c and d of ax^3 + bx^2 + cx + d
	std::vector<double> lp_r_p_function[4];		// [0-3] : a,b,c and d of ax^3 + bx^2 + cx + d
							// Convention ... for a range with knots : lp_r[i] ~ lp_r[i+1]
							// a,b,c,d : ..[0-3][i]

public:
		
	LonePair( const double frac_x, const double frac_y, const double frac_z, std::string species, std::string type, const Eigen::Vector3d* lattice_vector )
	: Atom(frac_x,frac_y,frac_z,species,type,lattice_vector)
	{
		this->lp_h_matrix.setZero();
		this->lp_h_matrix_tmp.setZero();
		for(int i=0;i<3;i++){ this->lp_h_matrix_derivatives[i].setZero(); }

		this->lp_gd.setZero();
		this->lp_gd_int.setZero();

		std::ifstream fp;
		std::string str;
		std::stringstream ss;
		double a,b,c,d;

		std::string knot_file, radial_s_file, radial_p_file;

		if( species == "Pb" )
		{	knot_file = RADIAL_Pb_KNOT;	radial_s_file = RADIAL_Pb_S;	radial_p_file = RADIAL_Pb_P;	}
		else if( species == "Sn" )
		{	knot_file = RADIAL_Sn_KNOT;	radial_s_file = RADIAL_Sn_S;	radial_p_file = RADIAL_Sn_P;	}
		else if( species == "Bi" )
		{	knot_file = RADIAL_Bi_KNOT;	radial_s_file = RADIAL_Bi_S;	radial_p_file = RADIAL_Bi_P;	}
		else
		{	std::cout << "No matching radial wavefunction data fount for lone pair species ...\n";
			std::exit(EXIT_FAILURE);
		}

		fp.open(knot_file);
		if( fp.is_open() )
		{
			std::getline(fp,str);
			ss << str;	ss >> this->lp_spline_knot_count;	// Get Number Of Knots
			for(int i=0;i<this->lp_spline_knot_count;i++)
			{	ss.clear(); ss.str("");
				std::getline(fp,str);
				ss << str;
				ss >> a;
				this->lp_r.push_back(a);
			}
			fp.close();
		}
		else
		{	std::cout << "Failed to read file : " << knot_file << std::endl;
			std::exit(EXIT_FAILURE);
		}

		fp.open(radial_s_file);
		if( fp.is_open() )
		{	std::getline(fp,str);
			for(int i=0;i<this->lp_spline_knot_count-1;i++)
			{	ss.clear(); ss.str("");
				std::getline(fp,str);
				ss << str;
				ss >> a >> b >> c >> d;
				
				this->lp_r_s_function[0].push_back(a);
				this->lp_r_s_function[1].push_back(b);
				this->lp_r_s_function[2].push_back(c);
				this->lp_r_s_function[3].push_back(d);
			}
			fp.close();
		}
		else
		{	std::cout << "Failed to read file : " << radial_s_file << std::endl;
			std::exit(EXIT_FAILURE);
		}

		fp.open(radial_p_file);
		if( fp.is_open() )
		{	std::getline(fp,str);
			for(int i=0;i<this->lp_spline_knot_count-1;i++)
			{	ss.clear(); ss.str("");
				std::getline(fp,str);
				ss << str;
				ss >> a >> b >> c >> d;
	
				this->lp_r_p_function[0].push_back(a);
				this->lp_r_p_function[1].push_back(b);
				this->lp_r_p_function[2].push_back(c);
				this->lp_r_p_function[3].push_back(d);
			}
			fp.close();
		}
		else
		{	std::cout << "Failed to read file : " << radial_p_file << std::endl;
			std::exit(EXIT_FAILURE);
		}

	}	// end Constructor

	// VirtualOverriding - must be in the same format
	virtual void SetFeature( const double core_charge, const double extra_charge=0, const double val1=0, const double val2=0 ) override
	{
		Atom::SetFeature(core_charge);	// Set LP core charge
		this->lp_charge = extra_charge;
		this->lp_lambda = val1;
	}

	virtual void ShowFrac() const override
	{	
		Eigen::Vector3d frac = Atom::GetCoreFrac();
		printf("%4.3s%8.4s%12.6lf%12.6lf%12.6lf%8.3lf",Atom::GetSpecies().c_str(),Atom::GetType().c_str(),frac(0),frac(1),frac(2),Atom::GetCoreCharge());
		printf("%12.3lf%10.3lf%12.9s%14.13s\n",this->lp_charge,this->lp_lambda,"LP(e)","/ Lambda(eV)");
	}

	virtual void ShowCart() const override
	{
		Eigen::Vector3d cart = Atom::GetCoreCart();
		printf("%4.3s%8.4s%12.6lf%12.6lf%12.6lf%8.3lf",Atom::GetSpecies().c_str(),Atom::GetType().c_str(),cart(0),cart(1),cart(2),Atom::GetCoreCharge());
		printf("%12.3lf%10.3lf%12.9s%14.13s\n",this->lp_charge,this->lp_lambda,"LP(e)","/ Lambda(eV)");
	}

	virtual void InitialiseDerivative() override
	{	
		Atom::InitialiseDerivative();	// derivative w.r.t. "lone" core
	}

	virtual void UpdateDerivativeInternal( const Eigen::Matrix3d& lattice_matrix )
	{
		Atom::UpdateDerivativeInternal( lattice_matrix );	// This only updates Atom private : Eigen::Vector3d cart_gd_int , i.e., internal derivative
	}

	virtual ~LonePair()
	{
		cnt++;
		std::cout << "LoneP Destructor ~" << cnt << std::endl;
	}
};

#endif
