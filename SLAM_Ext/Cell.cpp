#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <chrono>

#include <cstdlib>
#include <cmath>
#include <ctime>

#include "Cell.hpp"
#include "Manager.hpp"

std::string currentDateTime()
{
	std::time_t t = std::time(nullptr);
	std::tm* now = std::localtime(&t);

	char buffer[128];
	strftime(buffer, sizeof(buffer), "%m-%d-%Y %X", now);
	return buffer;
}

static int line_cnt = 0;

Cell::Cell( std::string input )
{
	Manager manager;

	// workload check variable

	this->NumberOfAtoms = 0;

	std::ifstream fp;

	fp.open(input);
	
	if( !fp.is_open() )
	{	std::cout << "Failed to read file ..." << std::endl;
		std::exit(EXIT_FAILURE);
	}
	else
	{
		int comment_loc;		// Comment locations 
		std::string str, tmp;		//
		std::stringstream ss;

		double frac[3];			// tmp load factional coordinate of atoms
		std::string species;		// species - e.g., Mg, O etc ...
		std::string type;		// type    - shell, core ?.. (type of describing a species)
		double charge;


		double a,b,c,al,be,ga;		// worksapce - real/reci lattice vector
		double v[3];			// workspace - real/reci lattice vector
	
		Eigen::Vector3d vv;		// workspace - vector operation

		for(std::string str; std::getline(fp,str); )
		{	line_cnt++;

			comment_loc = str.find('#');
			str = str.substr(0,comment_loc);
			tmp.clear();

			if( str.length() > 0 )
			{
				ss.clear(); ss.str("");
				ss << str;
				ss >> tmp;		

				if( !tmp.compare("cell") )
				{
					tmp.clear();

					std::getline(fp,str);	line_cnt++; 	// red one more line
					ss.clear(); ss.str("");			// clear 'ss' buffer

					ss << str;
					
					ss >> this->lattice_param(0);
					ss >> this->lattice_param(1);
					ss >> this->lattice_param(2);

					ss >> this->lattice_angle(0);	this->lattice_angle(0) = this->lattice_angle(0)/180.*M_PI;
					ss >> this->lattice_angle(1);	this->lattice_angle(1) = this->lattice_angle(1)/180.*M_PI;
					ss >> this->lattice_angle(2);	this->lattice_angle(2) = this->lattice_angle(2)/180.*M_PI;	// angle 'degree' -> 'radian'

					// PROCESSING REAL LATTICE VECTORS / http://gisaxs.com/index.php/Unit_cell
					v[0] = this->lattice_param(0);
					v[1] = 0.;
					v[2] = 0.;
					this->real_vector[0] << v[0],v[1],v[2];			// set real lattice vector a = (L1,0,0);

					v[0] = this->lattice_param(1)*std::cos(this->lattice_angle(2));
					v[1] = this->lattice_param(1)*std::sin(this->lattice_angle(2));
					v[2] = 0.;
					this->real_vector[1] << v[0],v[1],v[2];			// set real lattice vector b = (L2*Cos(gamma),L2*Sin(gamma),0);

					v[0] = this->lattice_param(2)*std::cos(this->lattice_angle(1));
					v[1] = this->lattice_param(2)*(std::cos(this->lattice_angle(0))-std::cos(this->lattice_angle(1))*std::cos(this->lattice_angle(2)))/std::sin(this->lattice_angle(2));
					v[2] = this->lattice_param(2)*std::sqrt(1.-pow(std::cos(this->lattice_angle(1)),2.)-pow(v[1]/this->lattice_param(2),2.));		
					this->real_vector[2] << v[0],v[1],v[2];			// set real lattice vector c

					// Setting LatticeMatrix - Eigen::Matrix3d
					//				a			b			c
					this->lattice_matrix << this->real_vector[0](0),this->real_vector[1](0),this->real_vector[2](0),	// {a1,b1,c1
								this->real_vector[0](1),this->real_vector[1](1),this->real_vector[2](1),	//  a2,b2,c2
								this->real_vector[0](2),this->real_vector[1](2),this->real_vector[2](2);	//  a3,b3,c3}

					this->lattice_matrix_transpose = this->lattice_matrix.transpose();
					this->lattice_matrix_inverse   = this->lattice_matrix.inverse();

					// Expected result : T_hkl = a h + b k + c l		
					//

					// PROCESSING RECIPROCAL LATTICE VECTORS
					vv = this->real_vector[1].cross(this->real_vector[2]);	// get b x c
					this->volume = this->real_vector[0].adjoint()*vv;	// get a.(b x c)

					// u vector
					vv = 2.*M_PI*this->real_vector[1].cross(this->real_vector[2])/this->volume;	// (b x c) / Vcell
					this->reci_vector[0] << vv(0), vv(1), vv(2);
					// v vector
					vv = 2.*M_PI*this->real_vector[2].cross(this->real_vector[0])/this->volume;	// (c x a) / Vcell
					this->reci_vector[1] << vv(0), vv(1), vv(2);
					// w vector
					vv = 2.*M_PI*this->real_vector[0].cross(this->real_vector[1])/this->volume;	// (a x b) / Vcell
					this->reci_vector[2] << vv(0), vv(1), vv(2);

					if( this->volume < 0. ) { this->volume = this->volume*-1.; }

					// Expected result : G_hkl = 2PI u h + 2PI v k + 2PI w l	// Note that 2.*M_PI factor is multiplied
					//
					//

				}// end if 'cell'

				if( !tmp.compare("lattice_vector") )
				{
					tmp.clear();
					std::getline(fp,str);	line_cnt++; 	// red one more line
					ss.clear(); ss.str("");			// clear 'ss' buffer
					ss << str;

					ss >> v[0] >> v[1] >> v[2];
					this->real_vector[0] << v[0],v[1],v[2];			// set real lattice vector a = (L1,0,0);
					this->lattice_param(0) = v[0];

					std::getline(fp,str);	line_cnt++; 	// red one more line
					ss.clear(); ss.str("");			// clear 'ss' buffer
					ss << str;
			
					ss >> v[0] >> v[1] >> v[2];
					this->real_vector[1] << v[0],v[1],v[2];			// set real lattice vector b = (L2*Cos(gamma),L2*Sin(gamma),0);
					this->lattice_param(1) = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
					this->lattice_angle(2) = acos(v[0]/this->lattice_param(1));

					std::getline(fp,str);	line_cnt++; 	// red one more line
					ss.clear(); ss.str("");			// clear 'ss' buffer
					ss << str;

					ss >> v[0] >> v[1] >> v[2];
					this->real_vector[2] << v[0],v[1],v[2];			// set real lattice vector c
					this->lattice_param(2) = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
					this->lattice_angle(1) = acos(v[0]/this->lattice_param(2));
					this->lattice_angle(0) = acos(v[1]/this->lattice_param(2)*sin(this->lattice_angle(2)) + cos(this->lattice_angle(1))*cos(this->lattice_angle(2)));

					// Setting LatticeMatrix - Eigen::Matrix3d
					//				a			b			c
					this->lattice_matrix << this->real_vector[0](0),this->real_vector[1](0),this->real_vector[2](0),	// {a1,b1,c1
								this->real_vector[0](1),this->real_vector[1](1),this->real_vector[2](1),	//  a2,b2,c2
								this->real_vector[0](2),this->real_vector[1](2),this->real_vector[2](2);	//  a3,b3,c3}
					this->lattice_matrix_transpose = this->lattice_matrix.transpose();
					this->lattice_matrix_inverse   = this->lattice_matrix.inverse();

					// PROCESSING RECIPROCAL LATTICE VECTORS
					vv = this->real_vector[1].cross(this->real_vector[2]);	// get b x c
					this->volume = this->real_vector[0].adjoint()*vv;	// get a.(b x c)

					// u vector
					vv = 2.*M_PI*this->real_vector[1].cross(this->real_vector[2])/this->volume;	// (b x c) / Vcell
					this->reci_vector[0] << vv(0), vv(1), vv(2);
					// v vector
					vv = 2.*M_PI*this->real_vector[2].cross(this->real_vector[0])/this->volume;	// (c x a) / Vcell
					this->reci_vector[1] << vv(0), vv(1), vv(2);
					// w vector
					vv = 2.*M_PI*this->real_vector[0].cross(this->real_vector[1])/this->volume;	// (a x b) / Vcell
					this->reci_vector[2] << vv(0), vv(1), vv(2);

					if( this->volume < 0. ) { this->volume = this->volume*-1.; }
				}

				if( !tmp.compare("atom") )
				{
					tmp.clear();
					ss >> frac[0] >> frac[1] >> frac[2] >> species >> type;

					if( !type.compare("core") )
					{	this->AtomList[NumberOfAtoms++] = new Atom(frac[0],frac[1],frac[2],species,type,this->real_vector);	// here 'real_vector' contains lattice vectors
					}
					else if( !type.compare("shel") )
					{
						double shel_frac[3];
						ss >> shel_frac[0] >> shel_frac[1] >> shel_frac[2];
						
						if( ss.fail() )
						{	
							this->AtomList[NumberOfAtoms++] = new Shell(frac[0],frac[1],frac[2],species,type,this->real_vector,frac[0],frac[1],frac[2]);
						}	// Exception 1 - if shell position is not specified
						else							
						{	
							this->AtomList[NumberOfAtoms++] = new Shell(frac[0],frac[1],frac[2],species,type,this->real_vector,shel_frac[0],shel_frac[1],shel_frac[2]);
							// Possibly add 'warning' if core-shell distance is too far
						}	// Exception 2 - if shell position is specified
					}
					else if( !type.compare("lone") )
					{	this->AtomList[NumberOfAtoms++] = new LonePair(frac[0],frac[1],frac[2],species,type,this->real_vector);	// here 'real_vector' contains lattice vectors
					}
					else	// TYPE CHECK
					{	std::cout << "error found in atom type ... ";
						std::cout << "specified type : " << "'" << type << "'" << " is not found in " << "'" << input << "'\n"; 
						std::cout << "list of allowed types : 'core', 'shel', and 'lone'\n";
						std::cout << "<SEE LINE> : " << line_cnt << std::endl;
						std::exit(EXIT_FAILURE);
					}
				}// end if 'atom'

				if( !tmp.compare("species") )
				{
					tmp.clear();
					ss >> species >> type >> charge;

					for(int i=0;i<this->NumberOfAtoms;i++)
					{	
						if( species == AtomList[i]->species && type == AtomList[i]->type )
						{
							if( !type.compare("core") )
							{	AtomList[i]->SetFeature(charge);
							}
							else if( !type.compare("shel") )
							{	
								double shel_charge, k2, k4;
								ss >> shel_charge >> k2 >> k4;
								AtomList[i]->SetFeature(charge,shel_charge,k2,k4);	// override function
							}
							else if( !type.compare("lone") )
							{	
								double lp_charge, lp_lambda;
								ss >> lp_charge >> lp_lambda;
								AtomList[i]->SetFeature(charge,lp_charge,lp_lambda);
							}
						}
					}

					// TYPE_CEHCK .. need implementation
					/*
						else	// TYPE CHECK
						{	std::cout << "no matching types found in species with atoms ... ";
							std::cout << "specified species " << "'" << species << "' ";
							std::cout << "or type " << "'" << type << "'" << " is not found in " << "'" << input << "' atom list\n"; 
							std::cout << "<SEE LINE> : " << line_cnt << std::endl;
							std::cout << "list of allowed types : 'core', 'shel', and 'lone'\n";
							std::exit(EXIT_FAILURE);
						}
					*/

				}// end if 'species'
			}//	end length
		}// end for getline()
	fp.close();
	}
	// Prepare Parameter ... for periodic sum
	manager.InitialisePeriodicSysParameter(*this);
}	// END 'Cell' class initilising


////	////	////	////

////	 Calculate Energy

////	////	////	////

void Cell::CalcCoulombEnergy()
{
	Manager manager;	// Managing class - interaction
	Eigen::Vector3d trans;
	Eigen::Vector3d delta_r, delta_rij;
	manager.InitialiseEnergy(*this);

	auto start = std::chrono::system_clock::now();

	for(int i=0;i<this->NumberOfAtoms;i++)
	{	for(int j=0;j<this->NumberOfAtoms;j++)
		{	for(int h = -this->h_max ; h <= this->h_max ; h++)
			{	for(int k = -this->k_max ; k <= this->k_max ; k++)
				{	for(int l = -this->l_max ; l <= this->l_max ; l++)
					{
						trans   = h*this->real_vector[0] + k*this->real_vector[1] + l*this->real_vector[2];	// T = h*a + k*b +l*c
						//// START REAL SPACE
						if( (trans.norm()) < this->rcut )
						{	
							if( h == 0 && k == 0 && l == 0 )
							{	
							    if( i != j )
							    {
								manager.CoulombMonoMonoReal(*this,i,j,trans);
							    }	// h=k=l=0 (central image) - excluding self interaction
							}
							else
							{
								manager.CoulombMonoMonoReal(*this,i,j,trans);
							}

							this->energy_real_sum_cnt++;
						}
						//// END REAL SPACE
					}// end l
				}// end k
			}//end h
		}// end j
	}// end i

	auto end = std::chrono::system_clock::now();

	this->energy_real_wtime = (end - start);

	start = std::chrono::system_clock::now();
  
	for(int i=0;i<this->NumberOfAtoms;i++)
	{	for(int j=0;j<this->NumberOfAtoms;j++)
		{	for(int h = -this->ih_max ; h <= this->ih_max ; h++)
			{	for(int k = -this->ik_max ; k <= this->ik_max ; k++)
				{	for(int l = -this->il_max ; l <= this->il_max ; l++)
					{
						trans   = h*this->reci_vector[0] + k*this->reci_vector[1] + l*this->reci_vector[2];	// G = h*2pi*u + k*2pi*v + l*2pi*w

						if( (trans.norm()) < this->gcut )
						{
							if( h == 0 && k == 0 && l == 0 )
							{
							    if( i == j )	// Self interaction
							    {	
								manager.CoulombMonoMonoSelf(*this,i,j,trans);
							    }
							}
							else
							{	
								manager.CoulombMonoMonoReci(*this,i,j,trans);
							}

							this->energy_reci_sum_cnt++;
						}
					}// end l
				}// end k
			}//end h
		}// end j
	}// end i

	// Wtime measure
	end = std::chrono::system_clock::now();
	this->energy_reci_wtime = (end - start);

	// Calculate Totl Energy
	this->mono_total_energy = this->mono_real_energy + this->mono_reci_energy + this->mono_reci_self_energy;

	return;
}

////	////	////	////

////	 Calculate Derivative

////	////	////	////

void Cell::CalcCoulombDerivative()
{
	Manager manager;
	Eigen::Vector3d trans;
	Eigen::Vector3d delta_rij, delta_r;
	manager.InitialiseDerivative(*this);

	auto start = std::chrono::system_clock::now();

	for(int i=0;i<this->NumberOfAtoms;i++)
	{	for(int j=0;j<this->NumberOfAtoms;j++)
		{	for(int h = -this->h_max ; h <= this->h_max ; h++)
			{	for(int k = -this->k_max ; k <= this->k_max ; k++)
				{	for(int l = -this->l_max ; l <= this->l_max ; l++)
					{
						//// START REAL SPACE
						trans   = h*this->real_vector[0] + k*this->real_vector[1] + l*this->real_vector[2];	// T = h*a + k*b +l*c

						if( (trans.norm()) < this->rcut )
						{	
							if( h == 0 && k == 0 && l == 0 )
							{	
							    if( i != j )
							    {	
								// Raw Geometric Derivatives
								manager.CoulombDerivativeReal(*this,i,j,trans);
								// Strain Derivatives Real Space Contribution
								manager.StrainDerivativeReal(*this,i,j,trans);
							    }	// h=k=l=0 (central image) - excluding self interaction
							}
							else
							{
								// Raw Geometric Derivatives
								manager.CoulombDerivativeReal(*this,i,j,trans);
								// Strain Derivatives Real Space Contribution
								manager.StrainDerivativeReal(*this,i,j,trans);
							}

							this->derivative_real_sum_cnt++;
						}
						//// END REAL SPACE
					}// end l
				}// end k
			}//end h
		}// end j
	}// end i
	auto end = std::chrono::system_clock::now();
	this->derivative_real_wtime = (end - start);

	start = std::chrono::system_clock::now();

	for(int i=0;i<this->NumberOfAtoms;i++)
	{	for(int j=0;j<this->NumberOfAtoms;j++)
		{	for(int h = -this->ih_max ; h <= this->ih_max ; h++)
			{	for(int k = -this->ik_max ; k <= this->ik_max ; k++)
				{	for(int l = -this->il_max ; l <= this->il_max ; l++)
					{
						//## START RECIPROCAL SPACE
						trans   = h*this->reci_vector[0] + k*this->reci_vector[1] + l*this->reci_vector[2];	// G = h*2pi*u + k*2pi*v + l*2pi*w
						
						if( (trans.norm()) < this->gcut )
						{
							if( h == 0 && k == 0 && l == 0 )
							{
								if( i == j )	// Self interaction
								{
									manager.CoulombDerivativeSelf(*this,i,j,trans);	// This calculation is not necessary for a system only with point charges
															// In the presence of ShellModel ions, there are non-zero derivatives by shell-core separation
									manager.StrainDerivativeSelf(*this,i,j,trans);
								}
							}
							else
							{	// Raw Geometric Derivative
								manager.CoulombDerivativeReci(*this,i,j,trans);
								// Strain Derivatives Reciprocal Space Contribution	// w.r.t 'r', 'V', 'G'
								manager.StrainDerivativeReci(*this,i,j,trans);
							}

							this->derivative_reci_sum_cnt++;
						}
					}// end l
				}// end k
			}//end h
		}// end j
	}// end i

	end = std::chrono::system_clock::now();
	this->derivative_reci_wtime = (end - start);

	// Convert raw geometric derivative -> internal geometric derivative
	for(int i=0;i<this->NumberOfAtoms;i++)
	{	
		AtomList[i]->UpdateDerivativeInternal(this->lattice_matrix_transpose);
		// see Note on ipad (Lattice Dynamics 25 July 2022)
	}
	// Symmeterise Strain Derivatives
	this->lattice_sd(1,0) = this->lattice_sd(0,1);
	this->lattice_sd(2,0) = this->lattice_sd(0,2);
	this->lattice_sd(2,1) = this->lattice_sd(1,2);

	return;
}

////	////	////	////

////	 LONEPAIR__ Calculate Energy

////	////	////	////

//#define LP_PBC_SUM_WTIME

void LonePair_begin_scf()
{
	using std::cout; using::std::endl;

	cout << "---------------------------------------------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "    Begin SCF Cycle" << endl;
	cout << endl;
	printf( "    SCF Tolerance (eV / per atom) : %10.4e\n",LONEPAIR_SCF_TOL);
	cout << endl;
	//cout << "---------------------------------------------------------------------------------------------------------" << endl;
	return;
}

void Cell::CalcLonePairCoulombEnergy()
{
	using std::cout; using std::endl;

	LonePair_begin_scf();
	bool is_first_scf = true;
	bool is_scf_done = false;

	Manager manager;
	Eigen::Vector3d trans;			// Translation Vector
	Eigen::Vector3d delta_r, delta_rij;	// Interatomic Distance

	int delimiter_cnt = 0;
	const int delimiter_max = 6;
	std::chrono::duration<double> pair_wtime;
	std::chrono::duration<double> real_wtime;
	std::chrono::duration<double> reci_wtime;

	/*	Initiate SCF Cycle
		1. Set Position Integral result (Real Space)
		2. Set Esp parameters + Diagonalisation (setting the ground-state)	*/
	manager.InitialiseSCF(*this);
	// This must be modified when optimiser comes together

	// SCF Loop 'n'
	for(int n=0;n<this->scf_iter_max;n++)
	{
		/*	Initialise(Zeros) for LPC_... LPLP_... temporal interaction values realted with LP entities
			Initialise(Zeros) Cell.lp_real_energy = Cell.lp_reci_energy = Cell.lp_reci_self_energy = Cell.lp_total_energy = 0;	*/
		manager.InitialiseLonePairCalculation_Energy(*this,is_first_scf);	// renew - lp_h_matrix_tmp with onsite terms (lp_lambda)

		// BLOCK OUTPUT
		cout << "---------------------------------------------------------------------------------------------------------" << endl;
		cout << endl;
		cout << "   SCF Cycle : " << std::noshowpos << n+1 << endl;
		cout << endl;
		cout << "---------------------------------------------------------------------------------------------------------" << endl;
		cout << endl;
		cout << "   Real space evaluation starts" << endl;
		cout << endl;
		// BLOCK OUTPUT

		auto real_sta = std::chrono::system_clock::now();

		for(int i=0;i<this->NumberOfAtoms;i++)
		{	for(int j=0;j<this->NumberOfAtoms;j++)
			{	
				auto pair_sta = std::chrono::system_clock::now();

				for(int h = -this->h_max ; h <= this->h_max ; h++)
				{	for(int k = -this->k_max ; k <= this->k_max ; k++)
					{	for(int l = -this->l_max ; l <= this->l_max ; l++)
						{
							trans   = h*this->real_vector[0] + k*this->real_vector[1] + l*this->real_vector[2];	// T = h*a + k*b +l*c
							//// START REAL SPACE
							if( (trans.norm()) < this->rcut )
							{	
								if( h == 0 && k == 0 && l == 0 )	// within the central sublattice
								{	
								    if( i != j )			// DO NOT directly compute self interaction !!
								    {
									manager.CoulombLonePairReal(*this,i,j,trans,is_first_scf);	// Sep 7 2022 Wed, for the moment calculating 'LonePair Core' involved interactions
																		// e.g., lp_core - core, lp_core - shel, lp_core - lp_core
									//printf("CENTRAL SUBLATTICE ----------- : %10.4d\t%10.4d\t%10.4d\n",h,k,l);
								    }	// h=k=l=0 (central image) - excluding self interaction
								}
								else					// case outside the central sublattice
								{
									//printf("non central : %10.4d\t%10.4d\t%10.4d\n",h,k,l);
									manager.CoulombLonePairReal(*this,i,j,trans,is_first_scf);
								}
							}
							//// END REAL SPACE
						}// end l
					}// end k
				}//end h

				auto pair_end = std::chrono::system_clock::now();
				pair_wtime = pair_end - pair_sta;
				// BLOCK OUTPUT
				printf("  (%d/%d)%6.3lfs",i+1,j+1,pair_wtime.count());
				// BLOCK OUTPUT
				if( delimiter_cnt == delimiter_max ){ delimiter_cnt = 0; cout << endl; }
				else { delimiter_cnt++; }

			}// end j
		}// end i

		delimiter_cnt = 0;

		auto real_end = std::chrono::system_clock::now();

		// BLOCK OUTPUT
		cout << endl;

		cout << endl;
		real_wtime = real_end - real_sta;
		cout << "   Real space wtime : " << real_wtime.count() << " s" << endl;
		cout << endl;
		cout << " *********************** " << endl;
		cout << endl;
		cout << "   Reciprocal space evaluation starts" << endl;
		cout << endl;
		// BLOCK OUTPUT

		auto reci_sta = std::chrono::system_clock::now();

		for(int i=0;i<this->NumberOfAtoms;i++)
		{	for(int j=0;j<this->NumberOfAtoms;j++)
			{
				auto pair_sta = std::chrono::system_clock::now();

				for(int h = -this->ih_max ; h <= this->ih_max ; h++)
				{	for(int k = -this->ik_max ; k <= this->ik_max ; k++)
					{	for(int l = -this->il_max ; l <= this->il_max ; l++)
						{
							trans   = h*this->reci_vector[0] + k*this->reci_vector[1] + l*this->reci_vector[2];	// G = h*2pi*u + k*2pi*v + l*2pi*w

							if( (trans.norm()) < this->gcut )
							{
								if( h == 0 && k == 0 && l == 0 )
								{
								    if( i == j )	// Self interaction
								    {	
									manager.CoulombLonePairSelf(*this,i,j,trans,is_first_scf);
								    }
								}
								else
								{	
									manager.CoulombLonePairReci(*this,i,j,trans,is_first_scf);
								}
							}
						}// end l
					}// end k
				}//end h

				auto pair_end = std::chrono::system_clock::now();
				pair_wtime = pair_end - pair_sta;
		
				// BLOCK OUTPUT
				printf("  (%d/%d)%6.3lfs",i+1,j+1,pair_wtime.count());
				// BLOCK OUTPUT
				
				if( delimiter_cnt == delimiter_max ){ delimiter_cnt = 0; cout << endl; }
				else { delimiter_cnt++; }

			}// end j
		}// end i
	
		// BLOCK OUTPUT
		delimiter_cnt = 0;

		auto reci_end = std::chrono::system_clock::now();
		cout << endl;

		cout << endl;
		reci_wtime = reci_end - reci_sta;
		cout << "   Reci space wtime : " << reci_wtime.count() << " s" << endl;
		cout << endl;
		cout << " *********************** " << endl;
		cout << endl;

		// Routine for Testify SCF Convergence
		// Diagonalisations : Determine GroundStates of LonePairs
		manager.GetLonePairGroundState(*this);

		cout << "   Lone Pair Information ('*' ground state)" << endl;
		cout << endl;
		for(int i=0;i<this->NumberOfAtoms;i++)
		{	if( this->AtomList[i]->type == "lone" )
			{
				cout << "    " << this->AtomList[i]->species << " (" << i+1 << ")" << endl;
				cout << endl;
				cout << "    Eigenvalues (eV) |  Eigenvectors [ s, px, py, pz ]" << endl;
				cout << endl;
				for(int j=0;j<4;j++)
				{
					if( j != static_cast<LonePair*>(this->AtomList[i])->GetGSIndex() )
					{	printf("   %s%14.9e  |\t",static_cast<LonePair*>(this->AtomList[i])->GetEval(j)>=0.?"+":"",static_cast<LonePair*>(this->AtomList[i])->GetEval(j));
						printf("%s%12.6e\t",static_cast<LonePair*>(this->AtomList[i])->GetEvec(0,j)>=0.?"+":"",static_cast<LonePair*>(this->AtomList[i])->GetEvec(0,j));
						printf("%s%12.6e\t",static_cast<LonePair*>(this->AtomList[i])->GetEvec(1,j)>=0.?"+":"",static_cast<LonePair*>(this->AtomList[i])->GetEvec(1,j));
						printf("%s%12.6e\t",static_cast<LonePair*>(this->AtomList[i])->GetEvec(2,j)>=0.?"+":"",static_cast<LonePair*>(this->AtomList[i])->GetEvec(2,j));
						printf("%s%12.6e\n",static_cast<LonePair*>(this->AtomList[i])->GetEvec(3,j)>=0.?"+":"",static_cast<LonePair*>(this->AtomList[i])->GetEvec(3,j));	}
					else
					{	printf(" * %s%14.9e  |\t",static_cast<LonePair*>(this->AtomList[i])->GetEval(j)>=0.?"+":"",static_cast<LonePair*>(this->AtomList[i])->GetEval(j));
						printf("%s%12.6e\t",static_cast<LonePair*>(this->AtomList[i])->GetEvec(0,j)>=0.?"+":"",static_cast<LonePair*>(this->AtomList[i])->GetEvec(0,j));
						printf("%s%12.6e\t",static_cast<LonePair*>(this->AtomList[i])->GetEvec(1,j)>=0.?"+":"",static_cast<LonePair*>(this->AtomList[i])->GetEvec(1,j));
						printf("%s%12.6e\t",static_cast<LonePair*>(this->AtomList[i])->GetEvec(2,j)>=0.?"+":"",static_cast<LonePair*>(this->AtomList[i])->GetEvec(2,j));
						printf("%s%12.6e\n",static_cast<LonePair*>(this->AtomList[i])->GetEvec(3,j)>=0.?"+":"",static_cast<LonePair*>(this->AtomList[i])->GetEvec(3,j));	}
				} 
				cout << endl;
			}
		}
		// BLOCK OUTPUT
		
		// only when it is the first SCF CYC - Calculating LonePair Core involved interactions
		if( is_first_scf == true )
		{	is_first_scf = false;	// false flag indicates, only calculate LonePair density
			this->mono_total_energy = this->mono_real_energy + this->mono_reci_energy + this->mono_reci_self_energy; // Update LonePair Core Contribution
		}

		// BLOCK OUTPUT
		cout << " *********************** " << endl;
		cout << endl;
		cout << "    Cycle Ends" << endl;
		cout << endl;
		cout << "    Cell Energy      (eV) : " << std::showpos << std::setprecision(16) << this->mono_total_energy + this->lp_total_energy << endl;
		cout << endl;
		cout << " *********************** " << endl;
		cout << endl;
		cout << "    (1) MM Energy Components" << endl;
		cout << endl;
		cout << "    Real              (eV) : " << std::showpos << std::setprecision(16) << this->mono_real_energy << endl;
		cout << "    Reciprocal        (eV) : " << std::showpos << std::setprecision(16) << this->mono_reci_energy + this->mono_reci_self_energy << endl;
		cout << "    MM Total Energy   (eV) : " << std::showpos << std::setprecision(16) << this->mono_total_energy << endl;
		cout << endl;
		cout << "    (2) LonePair Energy Components" << endl;
		cout << endl;
		cout << "    LP Real           (eV) : " << std::showpos << std::setprecision(16) << this->lp_real_energy << endl;
		cout << "    LP Reciprocal     (eV) : " << std::showpos << std::setprecision(16) << this->lp_reci_energy << endl;
		cout << "    LP Eigenvalue Sum (eV) : " << std::showpos << std::setprecision(16) << this->lp_eval_sum << "    (including self terms)" << endl;
		cout << endl;
		cout << "    LP Total Energy   (eV) : " << std::showpos << std::setprecision(16) << this->lp_total_energy << endl;
		cout << endl;
		// BLOCK OUTPUT

		if( manager.IsSCFDone( LONEPAIR_SCF_TOL/this->NumberOfAtoms ) ){ break; }
	}// end SCF

	// BLOCK OUTPUT
	cout << "---------------------------------------------------------------------------------------------------------" << endl;
	cout << endl;
	cout << "    Finalising SCF Process ..." << endl;
	cout << endl;
	cout << "    SCF Energy Profile (eV) [ Cell, LP Total, LP Eigenvalues, LP Real, LP Reciprocal ]" << endl;
	cout << endl;
	manager.PrintSCFProfile(*this);
	cout << endl;
	cout << "---------------------------------------------------------------------------------------------------------" << endl;
	// BLOCK OUTPUT

	return;
}


////	////	////	////

////	 LONEPAIR__ Calculate Derivative

////	////	////	////

void Cell::CalcLonePairCoulombDerivative()
{
	Manager manager;			// Managing Interactions + Calculations
	Eigen::Vector3d trans;			// Lattice Translation Vector
	Eigen::Vector3d delta_rij, delta_r;	// Interatomic Diatance

	manager.InitialiseLonePairCalculation_Derivatives(*this);

	auto start = std::chrono::system_clock::now();

	for(int i=0;i<this->NumberOfAtoms;i++)
	{	for(int j=0;j<this->NumberOfAtoms;j++)
		{	for(int h = -this->h_max ; h <= this->h_max ; h++)
			{	for(int k = -this->k_max ; k <= this->k_max ; k++)
				{	for(int l = -this->l_max ; l <= this->l_max ; l++)
					{
						//// START REAL SPACE
						trans   = h*this->real_vector[0] + k*this->real_vector[1] + l*this->real_vector[2];	// T = h*a + k*b +l*c

						if( (trans.norm()) < this->rcut )
						{	
							if( h == 0 && k == 0 && l == 0 )
							{	
							    if( i != j )
							    {	
								// Raw Geometric Derivatives
								manager.CoulombLonePairDerivativeReal(*this,i,j,trans);
								// Strain Derivatives Real Space Contribution
								manager.StrainLonePairDerivativeReal(*this,i,j,trans);
							    }	// h=k=l=0 (central image) - excluding self interaction
							}
							else
							{
								// Raw Geometric Derivatives
								manager.CoulombLonePairDerivativeReal(*this,i,j,trans);
								// Strain Derivatives Real Space Contribution
								manager.StrainLonePairDerivativeReal(*this,i,j,trans);
							}
						}
						//// END REAL SPACE
					}// end l
				}// end k
			}//end h
		}// end j
	}// end i
	auto end = std::chrono::system_clock::now();



	start = std::chrono::system_clock::now();

	for(int i=0;i<this->NumberOfAtoms;i++)
	{	for(int j=0;j<this->NumberOfAtoms;j++)
		{	for(int h = -this->ih_max ; h <= this->ih_max ; h++)
			{	for(int k = -this->ik_max ; k <= this->ik_max ; k++)
				{	for(int l = -this->il_max ; l <= this->il_max ; l++)
					{
						//## START RECIPROCAL SPACE
						trans   = h*this->reci_vector[0] + k*this->reci_vector[1] + l*this->reci_vector[2];	// G = h*2pi*u + k*2pi*v + l*2pi*w
						
						if( (trans.norm()) < this->gcut )
						{
							if( h == 0 && k == 0 && l == 0 )
							{
								if( i == j )	// Self interaction
								{
									manager.CoulombLonePairDerivativeSelf(*this,i,j,trans);	// This calculation is not necessary for a system only with point charges
															// In the presence of ShellModel ions, there are non-zero derivatives by shell-core separation
									manager.StrainLonePairDerivativeSelf(*this,i,j,trans);
								}
							}
							else
							{	// Raw Geometric Derivative
								manager.CoulombLonePairDerivativeReci(*this,i,j,trans);
								// Strain Derivatives Reciprocal Space Contribution	// w.r.t 'r', 'V', 'G'
								manager.StrainLonePairDerivativeReci(*this,i,j,trans);
							}
						}
					}// end l
				}// end k
			}//end h
		}// end j
	}// end i

	end = std::chrono::system_clock::now();

	// Convert raw geometric derivative -> internal geometric derivative
	for(int i=0;i<this->NumberOfAtoms;i++)
	{	
		AtomList[i]->UpdateDerivativeInternal(this->lattice_matrix_transpose);
		// see Note on ipad (Lattice Dynamics 25 July 2022)
	}
	// Symmeterise Strain Derivatives
	this->lattice_sd(1,0) = this->lattice_sd(0,1);
	this->lattice_sd(2,0) = this->lattice_sd(0,2);
	this->lattice_sd(2,1) = this->lattice_sd(1,2);

	return;
}

////	////	////	////

//// Calculate Lattice Gradient

////	////	////	////

void Cell::CalcLatticeDerivative() 	// Independe to LonePair Derivatives ... as long as the Lattice Strain Derivatives are measured
{
	this->lattice_derivative.setZero();

	// LatticeVector a = ( a1 , a2 , a3 ) L = 1,2,3
	for(int i=0;i<3;i++)	// loop for 'u'
	{	
		// LatticeVector a = ( a1 , a2 , a3 ) L = 1,2,3
		this->lattice_derivative(0,0) += this->lattice_sd(0,i)*this->lattice_matrix_inverse(0,i);
		this->lattice_derivative(1,0) += this->lattice_sd(1,i)*this->lattice_matrix_inverse(0,i);
		this->lattice_derivative(2,0) += this->lattice_sd(2,i)*this->lattice_matrix_inverse(0,i);
		//			 L n			  L u				    n u			    

		// LatticeVector b = ( b1 , b2 , b3 ) L = 1,2,3
		this->lattice_derivative(0,1) += this->lattice_sd(0,i)*this->lattice_matrix_inverse(1,i);
		this->lattice_derivative(1,1) += this->lattice_sd(1,i)*this->lattice_matrix_inverse(1,i);
		this->lattice_derivative(2,1) += this->lattice_sd(2,i)*this->lattice_matrix_inverse(1,i);

		// LatticeVector c = ( c1 , c2 , c3 ) L = 1,2,3
		this->lattice_derivative(0,2) += this->lattice_sd(0,i)*this->lattice_matrix_inverse(2,i);
		this->lattice_derivative(1,2) += this->lattice_sd(1,i)*this->lattice_matrix_inverse(2,i);
		this->lattice_derivative(2,2) += this->lattice_sd(2,i)*this->lattice_matrix_inverse(2,i);
	}
}




////	////	////	////

////	Cell Destructor

////	////	////	////


Cell::~Cell()
{
	for(auto i=0;i<this->NumberOfAtoms;i++)
	{	
		delete AtomList[i];
	}
}

////	////	////	////

////	 OUTPUT CHECK

////	////	////	////

void Cell::ShowBasicCellInfo() const
{
	using std::cout;
	using std::endl;

	cout << endl;
	cout << "*********************************************************************************************************\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << " Basic Cell Information\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << endl;
	printf("     Lattice Parameters (Angs) : %12.6lf\t%12.6lf\t%12.6lf\n",this->lattice_param(0),this->lattice_param(1),this->lattice_param(2));
	printf("     Lattice Angles     (Degs) : %12.6lf\t%12.6lf\t%12.6lf\n",180/M_PI*this->lattice_angle(0),180/M_PI*this->lattice_angle(1),180/M_PI*this->lattice_angle(2));
	cout << endl;
	printf("     Cell Volume      (Angs^3) : %12.6lf\n",this->volume);
	cout << endl;
	cout <<"                   x           y          z           (Transposed Lattice Matrix)" << endl;
	printf("       a | %12.6lf%12.6lf%12.6lf\n",real_vector[0](0),real_vector[0](1),real_vector[0](2));
	printf("       b | %12.6lf%12.6lf%12.6lf\n",real_vector[1](0),real_vector[1](1),real_vector[1](2));
	printf("       c | %12.6lf%12.6lf%12.6lf\n",real_vector[2](0),real_vector[2](1),real_vector[2](2));
	cout << endl;
	cout <<"                   u           v          w\n";
	printf("       u | %12.6lf%12.6lf%12.6lf\n",reci_vector[0](0)/(2.*M_PI),reci_vector[0](1)/(2.*M_PI),reci_vector[0](2)/(2.*M_PI));
	printf("       v | %12.6lf%12.6lf%12.6lf\n",reci_vector[1](0)/(2.*M_PI),reci_vector[1](1)/(2.*M_PI),reci_vector[1](2)/(2.*M_PI));
	printf("       w | %12.6lf%12.6lf%12.6lf\n",reci_vector[2](0)/(2.*M_PI),reci_vector[2](1)/(2.*M_PI),reci_vector[2](2)/(2.*M_PI));
	cout << endl;
	printf(" 2pi * u | %12.6lf%12.6lf%12.6lf\n",reci_vector[0](0),reci_vector[0](1),reci_vector[0](2));
	printf(" 2pi * v | %12.6lf%12.6lf%12.6lf\n",reci_vector[1](0),reci_vector[1](1),reci_vector[1](2));
	printf(" 2pi * w | %12.6lf%12.6lf%12.6lf\n",reci_vector[2](0),reci_vector[2](1),reci_vector[2](2));
	cout << endl;
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << "    Atom fractional coordinate\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << "  species           x           y           z      core_charge\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	for(auto i=0;i<NumberOfAtoms;i++)
	{	AtomList[i]->ShowFrac();
	}
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << "    Atom Cartesian coordinate\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << "  species           x           y           z      core_charge\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	for(auto i=0;i<NumberOfAtoms;i++)
	{	AtomList[i]->ShowCart();
	}
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << "*********************************************************************************************************\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << endl;
	cout << "    Lattice Summation Method (Ewald) / Accuracy Factor : " << -log10(this->accuracy) << endl;
	cout << endl;
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << endl;
	cout << "    Sigma  (Angs)   : " << std::setprecision(6) << this->sigma  << endl;
	cout << "    Weight          : " << std::setprecision(6) << this->weight << endl;
	cout << "    Rcut   (Angs)   : " << std::setprecision(6) << this->rcut   << endl;
	cout << "    Gcut   (Angs-1) : " << std::setprecision(6) << this->gcut   << endl;
	cout << "    hkl    (real)   : " << std::setprecision(6) << this->h_max << " " << this->k_max << " " << this->l_max << endl;
	cout << "    hkl    (reci)   : " << std::setprecision(6) << this->ih_max << " " << this->ik_max << " " << this->il_max << endl;
	cout << endl;
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << "*********************************************************************************************************\n";
	//cout << endl;
}


void Cell::ShowEnergyDerivative() const
{
	using std::cout;
	using std::endl;

	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << endl;
	cout << "    Total Component of energy :\n";
	cout << endl;
	//cout << "    Reci Self         (eV) : " << std::setprecision(16) << this->mono_reci_self_energy << endl;
	//cout << "    Reci              (eV) : " << std::setprecision(16) << this->mono_reci_energy << endl;
	cout << "    Real              (eV) : " << std::setprecision(16) << this->mono_real_energy + this->lp_real_energy << endl;
	//cout << "    Reciprocal + Self (eV) : " << std::setprecision(16) << this->mono_reci_self_energy + this->mono_reci_energy << endl;
	cout << "    Reciprocal        (eV) : " << std::setprecision(16) << this->mono_reci_self_energy + this->mono_reci_energy + this->lp_reci_energy << endl;
	cout << "    LonePair Evals    (eV) : " << std::setprecision(16) << this->lp_eval_sum << endl;
	cout << "    Total             (eV) : " << std::setprecision(16) << this->mono_total_energy + this->lp_total_energy << std::endl;
	cout << endl;
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << "    Geometric Derivatives            (eV/Angs)\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	for(int i=0;i<this->NumberOfAtoms;i++)
	{	printf("%6.4s\t%12.6lf\t%12.6lf\t%12.6lf\n",this->AtomList[i]->species.c_str(),this->AtomList[i]->cart_gd(0),this->AtomList[i]->cart_gd(1),this->AtomList[i]->cart_gd(2));
	}
	for(int i=0;i<this->NumberOfAtoms;i++)
	{	if( !AtomList[i]->type.compare("shel") )
		{
		printf("%6.4s\t%12.6lf\t%12.6lf\t%12.6lf\n",this->AtomList[i]->species.c_str(),
				static_cast<Shell*>(this->AtomList[i])->shel_cart_gd(0),
				static_cast<Shell*>(this->AtomList[i])->shel_cart_gd(1),
				static_cast<Shell*>(this->AtomList[i])->shel_cart_gd(2));
		}
	}
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << "    Internal Geometric Derivatives   (eV)\n"; // Lattice Matrix * Grad(i) E\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	for(int i=0;i<this->NumberOfAtoms;i++)
	{	printf("%6.4s\t%12.6lf\t%12.6lf\t%12.6lf\n",this->AtomList[i]->species.c_str(),this->AtomList[i]->cart_gd_int(0),this->AtomList[i]->cart_gd_int(1),this->AtomList[i]->cart_gd_int(2));
	}
	for(int i=0;i<this->NumberOfAtoms;i++)
	{	if( !AtomList[i]->type.compare("shel") )
		{
		printf("%6.4s\t%12.6lf\t%12.6lf\t%12.6lf\n",this->AtomList[i]->species.c_str(),
				static_cast<Shell*>(this->AtomList[i])->shel_cart_gd_int(0),
				static_cast<Shell*>(this->AtomList[i])->shel_cart_gd_int(1),
				static_cast<Shell*>(this->AtomList[i])->shel_cart_gd_int(2));
		}
	}
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << "    Strain Derivatives               (eV)\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	printf("%12.6s%12.6s%12.6s","e1(xx)","e6(xy)","e5(xz)"); printf("    %12.6lf\t%12.6lf\t%12.6lf\n",this->lattice_sd(0,0),this->lattice_sd(0,1),this->lattice_sd(0,2));
	printf("%12.6s%12.6s%12.6s","e6(yx)","e2(yy)","e4(yz)"); printf("    %12.6lf\t%12.6lf\t%12.6lf\n",this->lattice_sd(1,0),this->lattice_sd(1,1),this->lattice_sd(1,2));
	printf("%12.6s%12.6s%12.6s","e5(zx)","e4(zy)","e3(zz)"); printf("    %12.6lf\t%12.6lf\t%12.6lf\n",this->lattice_sd(2,0),this->lattice_sd(2,1),this->lattice_sd(2,2));
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << "    Lattice Derivatives               (eV/Angs)\n";
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << endl;
	cout << "             a           b          c " << endl;
	printf( "%5.3s%12.6lf%12.6lf%12.6lf\n","x",this->lattice_derivative(0,0),this->lattice_derivative(0,1),this->lattice_derivative(0,2));
	printf( "%5.3s%12.6lf%12.6lf%12.6lf\n","y",this->lattice_derivative(1,0),this->lattice_derivative(1,1),this->lattice_derivative(1,2));
	printf( "%5.3s%12.6lf%12.6lf%12.6lf\n","z",this->lattice_derivative(2,0),this->lattice_derivative(2,1),this->lattice_derivative(2,2));
	cout << endl;
	cout << "    Lattice Matrix (reference) " << endl;
	cout << endl;
	cout << "             a           b          c " << endl;
	printf( "%5.3s%12.6lf%12.6lf%12.6lf\n","x",this->lattice_matrix(0,0),this->lattice_matrix(0,1),this->lattice_matrix(0,2));
	printf( "%5.3s%12.6lf%12.6lf%12.6lf\n","y",this->lattice_matrix(1,0),this->lattice_matrix(1,1),this->lattice_matrix(1,2));
	printf( "%5.3s%12.6lf%12.6lf%12.6lf\n","z",this->lattice_matrix(2,0),this->lattice_matrix(2,1),this->lattice_matrix(2,2));
	cout << endl;
	cout << "    Cell energy : "; printf("%12.12lf",this->mono_total_energy + this->lp_total_energy); cout << "  (eV)" << endl;
	cout << endl;
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << endl;
/*
	cout << endl;
	cout << " Calculation Wtime ( Energy )\n";
	cout << endl;
	cout << " Real           : " << std::setprecision(6) << this->energy_real_wtime.count() << " s\n";
	cout << " Reciprocal Sum : " << std::setprecision(6) << this->energy_reci_wtime.count() << " s\n";
	cout << endl;
	cout << " Calculation Wtime ( Derivatives )\n";
	cout << endl;
	cout << " Real           : " << std::setprecision(6) << this->derivative_real_wtime.count() << " s\n";
	cout << " Reciprocal Sum : " << std::setprecision(6) << this->derivative_reci_wtime.count() << " s\n";
	cout << endl;
	cout << " Number of Operations for Real/Reciprocal Summations\n";
	cout << endl;
	cout << " Energy     Real/Reciprocal : " << this->energy_real_sum_cnt << " / " << this->energy_reci_sum_cnt << endl;
	cout << " Derivative Real/Reciprocal : " << this->derivative_real_sum_cnt << " / " << this->derivative_reci_sum_cnt << endl;
	cout << endl;
	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << endl;
*/
/*
	cout << this->AtomList[0]->GetCoreCharge() << endl;
	if( this->AtomList[0]->type == "shel" )
		cout << static_cast<Shell*>(this->AtomList[0])->GetShellCharge() << endl;
*/
}

void Cell::Finalise() const
{
using std::cout;
using std::endl;

	cout << "---------------------------------------------------------------------------------------------------------\n";
	cout << " Job Finished at                  : " << currentDateTime() << endl;
	cout << "---------------------------------------------------------------------------------------------------------\n";

}
