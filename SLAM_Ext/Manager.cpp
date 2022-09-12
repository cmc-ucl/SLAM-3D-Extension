#include "Manager.hpp"

//#define MANAGER_DEBUG_lp_add_2
//#define MD2
#define MD3

// Implement Integrators

const Eigen::Matrix4d& Manager::set_h_matrix_real_pc( /* IN/RES OUT */ LonePair* lp, const Eigen::Vector3d& R, const double sig, const int lp_i, const int pc_i )
{
	this->man_lp_matrix_h.GetTransformationMatrix(R);	// get Transformation matrix ... saved : Eigen::Matrix4d this->man_lp_matrix_h.transform_matrix;
	Eigen::Matrix4d h_tmp;					// calculating local h_matrix WS
	h_tmp.setZero();					// initialise

	// evalulation block
	h_tmp(0,0) = this->man_lp_matrix_h.real_ss_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());
	h_tmp(0,3) = this->man_lp_matrix_h.real_sz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());
	h_tmp(3,0) = h_tmp(0,3);
	h_tmp(1,1) = this->man_lp_matrix_h.real_xx_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());
	h_tmp(2,2) = h_tmp(1,1);
	h_tmp(3,3) = this->man_lp_matrix_h.real_zz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());

	this->real_lp_h_pc[lp_i][pc_i] = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp * this->man_lp_matrix_h.transform_matrix;	// inverse transformation

	return this->real_lp_h_pc[lp_i][pc_i];
/*
	// Test 
	//cout << h_tmp << endl;
	for(int m=0;m<4;m++)
	{	for(int v=0;v<4;v++)
		{	for(int h=0;h<4;h++)
			{	for(int l=0;l<4;l++)
				{	//lp->lp_h_matrix_tmp(m,v) += this->man_lp_matrix_h.transform_matrix(h,m)*this->man_lp_matrix_h.transform_matrix(l,v)*h_tmp(h,l);
					this->real_lp_h_pc[lp_i][pc_i](m,v) +=  this->man_lp_matrix_h.transform_matrix(h,m)*this->man_lp_matrix_h.transform_matrix(l,v)*h_tmp(h,l);
				}
			}
		}
	}

	cout << "Final Res" << endl;
	cout << this->real_lp_h_pc[lp_i][pc_i] << endl;

	// Inverse Transformation - H_global = P_transpose * H_local * P
	cout << "Inverse Transformation ... Using eigen feature" << endl;
	Eigen::Matrix4d res = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp * this->man_lp_matrix_h.transform_matrix;
	cout << res << endl;

	// Direct Transformation - H_local = P * H_local * P_transpose
	cout << "Direct Transformation ... Using eigen feature" << endl;
	Eigen::Matrix4d res2 = this->man_lp_matrix_h.transform_matrix * res  * this->man_lp_matrix_h.transform_matrix.transpose();
	cout << res2 << endl;
*/
}

void Manager::set_h_matrix_real_pc_derivative( LonePair* lp, const Eigen::Vector3d& R, const double sig, const int lp_i, const int pc_i )
{
	this->man_lp_matrix_h.GetTransformationMatrix(R);	// get Transformation matrix ... saved : Eigen::Matrix4d this->man_lp_matrix_h.transform_matrix;
	Eigen::Matrix4d h_tmp_x, h_tmp_y, h_tmp_z;		// calculating local h_matrix WS
	Eigen::Matrix4d h_tmp_x_ws,h_tmp_y_ws,h_tmp_z_ws;
	Eigen::Vector3d v_loc, v_glo;

	h_tmp_x.setZero();					// initialise
	h_tmp_y.setZero();					// initialise
	h_tmp_z.setZero();					// initialise
	h_tmp_x_ws.setZero();					// initialise
	h_tmp_y_ws.setZero();					// initialise
	h_tmp_z_ws.setZero();					// initialise
	
	// evalulation block
	h_tmp_x(0,1) = this->man_lp_matrix_h.real_sx_grad_x_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());
	h_tmp_x(1,0) = h_tmp_x(0,1);
	h_tmp_x(1,3) = this->man_lp_matrix_h.real_xz_grad_x_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());
	h_tmp_x(3,1) = h_tmp_x(1,3);

	h_tmp_y(0,2) = h_tmp_y(2,0) = h_tmp_x(0,1);	// y-sy = x-sx
	h_tmp_y(2,3) = h_tmp_y(3,2) = h_tmp_x(1,3);	// y-yz = x-sz
	
	h_tmp_z(0,0) = this->man_lp_matrix_h.real_ss_grad_z_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());
	h_tmp_z(0,3) = this->man_lp_matrix_h.real_sz_grad_z_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());
	h_tmp_z(3,0) = h_tmp_z(0,3);
	h_tmp_z(1,1) = this->man_lp_matrix_h.real_xx_grad_z_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());
	h_tmp_z(2,2) = h_tmp_z(1,1);
	h_tmp_z(3,3) = this->man_lp_matrix_h.real_zz_grad_z_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());

	// note : h_tmp_*_ws are in the local reference frame, their x'/y'/z' element (local) has to be inversed to x/y/z (global) 
	h_tmp_x_ws = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp_x * this->man_lp_matrix_h.transform_matrix;
	h_tmp_y_ws = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp_y * this->man_lp_matrix_h.transform_matrix;
	h_tmp_z_ws = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp_z * this->man_lp_matrix_h.transform_matrix;

	for(int i=0;i<4;i++)
	{	for(int j=0;j<4;j++)
		{
			v_loc << h_tmp_x_ws(i,j), h_tmp_y_ws(i,j), h_tmp_z_ws(i,j);
			v_glo = this->man_lp_matrix_h.transform_matrix_shorthand.transpose() * v_loc;

			this->real_lp_h_pc_x[lp_i][pc_i](i,j) =  v_glo(0);
			this->real_lp_h_pc_y[lp_i][pc_i](i,j) =  v_glo(1);
			this->real_lp_h_pc_z[lp_i][pc_i](i,j) =  v_glo(2);
		}
	}
}
/*
void Manager::set_h_matrix_real_pc_derivative_y( LonePair* lp, const Eigen::Vector3d& R, const double sig, const int lp_i, const int pc_i )
{
	this->man_lp_matrix_h.GetTransformationMatrix(R);	// get Transformation matrix ... saved : Eigen::Matrix4d this->man_lp_matrix_h.transform_matrix;
	Eigen::Matrix4d h_tmp;					// calculating local h_matrix WS
	h_tmp.setZero();					// initialise
	
	// evalulation block
	h_tmp(0,2) = this->man_lp_matrix_h.real_sx_grad_x_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());
	h_tmp(2,0) = h_tmp(0,2);
	h_tmp(2,3) = this->man_lp_matrix_h.real_xz_grad_x_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());
	h_tmp(3,2) = h_tmp(2,3);
	
	this->real_lp_h_pc_y[lp_i][pc_i] = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp * this->man_lp_matrix_h.transform_matrix;	// inverse transformation
}

void Manager::set_h_matrix_real_pc_derivative_z( LonePair* lp, const Eigen::Vector3d& R, const double sig, const int lp_i, const int pc_i )
{
	this->man_lp_matrix_h.GetTransformationMatrix(R);	// get Transformation matrix ... saved : Eigen::Matrix4d this->man_lp_matrix_h.transform_matrix;
	Eigen::Matrix4d h_tmp;					// calculating local h_matrix WS
	h_tmp.setZero();					// initialise
	
	// evalulation block
	h_tmp(0,0) = this->man_lp_matrix_h.real_ss_grad_z_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());
	h_tmp(0,3) = this->man_lp_matrix_h.real_sz_grad_z_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());
	h_tmp(3,0) = h_tmp(0,3);
	h_tmp(1,1) = this->man_lp_matrix_h.real_xx_grad_z_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());
	h_tmp(2,2) = h_tmp(1,1);
	h_tmp(3,3) = this->man_lp_matrix_h.real_zz_grad_z_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());
	
	this->real_lp_h_pc_z[lp_i][pc_i] = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp * this->man_lp_matrix_h.transform_matrix;	// inverse transformation
}
*/

void Manager::InitialiseEnergy( Cell& C )
{
#ifdef MD3
	using std::cout, std::endl;
	LonePair* lp = nullptr;
	int ind;

	for(int i=0;i<C.NumberOfAtoms;i++)
	{
		cout << i << "\t" << C.AtomList[i]->type << endl;

		if( C.AtomList[i]->type == "lone" )
		{
			cout << "LonePairFound" << endl;
			ind = i;
			lp = static_cast<LonePair*>(C.AtomList[i]);
		}
	}
	double delta = 0.005;
	double sig =2.3231241;
	Eigen::Vector3d v;
	v << 1,2,-3;
	Manager::set_h_matrix_real_pc( lp, v, sig, ind, 0 );
	//cout << this->real_lp_h_pc[ind][0] << endl;
	for(int i=0;i<4;i++)
	{	for(int j=0;j<4;j++)
			printf("%20.12lf\t",this->real_lp_h_pc[ind][0](i,j));
		cout << endl;
	}

	Eigen::Matrix4d fdm_x[2];
	Eigen::Matrix4d res_x;
	res_x.setZero();
	fdm_x[0].setZero();
	fdm_x[1].setZero();
	Eigen::Matrix4d fdm_y[2];
	Eigen::Matrix4d res_y;
	res_y.setZero();
	fdm_y[0].setZero();
	fdm_y[1].setZero();
	Eigen::Matrix4d fdm_z[2];
	Eigen::Matrix4d res_z;
	res_z.setZero();
	fdm_z[0].setZero();
	fdm_z[1].setZero();

	cout << " fdm real_h_pc_x : " << delta << endl;
	v(0) = v(0) + delta;
	fdm_x[0] = Manager::set_h_matrix_real_pc( lp, v, sig, ind, 0 );
	v(0) = v(0) - delta;
	v(0) = v(0) - delta;
	fdm_x[1] = Manager::set_h_matrix_real_pc( lp, v, sig, ind, 0 );
	v(0) = v(0) + delta;
	res_x = (fdm_x[0] - fdm_x[1])/2./delta;
	for(int i=0;i<4;i++)
	{	for(int j=0;j<4;j++)
			printf("%20.12lf\t",res_x(i,j));
		cout << endl;
	}
	cout << " fdm real_h_pc_y : " << delta << endl;
	v(1) = v(1) + delta;
	fdm_y[0] = Manager::set_h_matrix_real_pc( lp, v, sig, ind, 0 );
	v(1) = v(1) - delta;
	v(1) = v(1) - delta;
	fdm_y[1] = Manager::set_h_matrix_real_pc( lp, v, sig, ind, 0 );
	v(1) = v(1) + delta;
	res_y = (fdm_y[0] - fdm_y[1])/2./delta;
	for(int i=0;i<4;i++)
	{	for(int j=0;j<4;j++)
			printf("%20.12lf\t",res_y(i,j));
		cout << endl;
	}
	cout << " fdm real_h_pc_z : " << delta << endl;
	v(2) = v(2) + delta;
	fdm_z[0] = Manager::set_h_matrix_real_pc( lp, v, sig, ind, 0 );
	v(2) = v(2) - delta;
	v(2) = v(2) - delta;
	fdm_z[1] = Manager::set_h_matrix_real_pc( lp, v, sig, ind, 0 );
	v(2) = v(2) + delta;
	res_z = (fdm_z[0] - fdm_z[1])/2./delta;
	for(int i=0;i<4;i++)
	{	for(int j=0;j<4;j++)
			printf("%20.12lf\t",res_z(i,j));
		cout << endl;
	}

	cout << "End of H Matrix Test" << endl;
	cout << endl;
	
	//Manager::set_h_matrix_real_pc_derivative_x( lp, v, sig, ind ,0 );
	//Manager::set_h_matrix_real_pc_derivative_y( lp, v, sig, ind ,0 );
	//Manager::set_h_matrix_real_pc_derivative_z( lp, v, sig, ind ,0 );
	Manager::set_h_matrix_real_pc_derivative( lp, v, sig, ind ,0 );
	
	cout << "grad h_pc_x" << endl;
	for(int i=0;i<4;i++)
	{	for(int j=0;j<4;j++)
			printf("%20.12lf\t",this->real_lp_h_pc_x[ind][0](i,j));
		cout << endl;
	}
	cout << "grad h_pc_y" << endl;
	for(int i=0;i<4;i++)
	{	for(int j=0;j<4;j++)
			printf("%20.12lf\t",this->real_lp_h_pc_y[ind][0](i,j));
		cout << endl;
	}
	cout << "grad h_pc_z" << endl;
	for(int i=0;i<4;i++)
	{	for(int j=0;j<4;j++)
			printf("%20.12lf\t",this->real_lp_h_pc_z[ind][0](i,j));
		cout << endl;
	}
	// FDM TEST 

	exit(1);
#endif
	
#ifdef MANAGER_DEBUG_lp_add_2
	// Test Section Sep 07 2022
	using std::cout, std::endl;

	Eigen::Vector3d v;
	v << 1,2,-3;
	cout << "Vector Size : " << v.norm() << endl;

	//cout << this->man_lp_matrix_h.transform_matrix << endl;
	this->man_lp_matrix_h.GetTransformationMatrix(v);
	printf("vector : %12.6lf\t%12.6lf\t%12.6lf\n",v(0),v(1),v(2));

	cout << endl;
	cout << "Transform matrix 4x4 version" << endl;
	cout << this->man_lp_matrix_h.transform_matrix << endl;
	cout << endl;
	cout << "Transform matrix 3x3 version" << endl;
	cout << this->man_lp_matrix_h.transform_matrix_shorthand << endl;
	
	cout << "Test Transformation\n";

	Eigen::Vector3d v2;
	/*
	v2.setZero();
	for(int i=0;i<3;i++)
	{	for(int j=0;j<3;j++)
		{	v2(i) += this->man_lp_matrix_h.transform_matrix(i+1,j+1)*v(j);
		}
	}	// v2 is transformed 
	*/
	v2 = this->man_lp_matrix_h.transform_matrix_shorthand * v;
	printf("vector : %12.6lf\t%12.6lf\t%12.6lf\n",v2(0),v2(1),v2(2));

	cout << endl;
	cout << "Inverse Test" << endl;
	Eigen::Vector3d v3;
	v3.setZero();
	/*
	for(int i=0;i<3;i++)
	{	for(int j=0;j<3;j++)
		{
			v3(i) += this->man_lp_matrix_h.transform_matrix(j+1,i+1)*v2(j);
		}
	}
	*/
	v3 = this->man_lp_matrix_h.transform_matrix_shorthand.transpose() * v2;
	printf("vector : %12.6lf\t%12.6lf\t%12.6lf\n",v3(0),v3(1),v3(2));

	cout << "reflected" << endl;
	for(int i=0;i<C.NumberOfAtoms;i++)
	{
		if( C.AtomList[i]->type == "lone" )
		{
			cout << "LonePairFound" << endl;
			LonePair* lp = static_cast<LonePair*>(C.AtomList[i]);
			this->man_lp_matrix_h.NIntegral_test_real( lp->lp_r, lp->lp_r_s_function, lp->lp_r_p_function );

		}
	}
	exit(1);
	// ForceQuit - 1
	// Test Section End
#endif

#ifdef MD2
	double sig  = 2.49475;
	//double dist = 1.305;
	double dist = 1.355;
	double res_ss, res_sz, res_xx, res_zz;
	using std::cout, std::endl;

	LonePair* lp;

	for(int i=0;i<C.NumberOfAtoms;i++)
	{
		if( C.AtomList[i]->type == "lone" )
		{
			cout << "LonePairFound" << endl;
			lp = static_cast<LonePair*>(C.AtomList[i]);
			res_ss = this->man_lp_matrix_h.real_ss_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,dist);
			res_sz = this->man_lp_matrix_h.real_sz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,dist);
			res_xx = this->man_lp_matrix_h.real_xx_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,dist);
			res_zz = this->man_lp_matrix_h.real_zz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,dist);
		}
	}
	printf("%10.6lf%20.10lf\t%20.10lf\t%20.10lf\t%20.10lf\n",dist,res_ss,res_sz,res_xx,res_zz);

	
	double delta[6] = {0.1,0.05,0.025,0.0125,0.005,0.00025};
	double ssf, ssb;

	cout << endl;
	for(int i=0;i<6;i++)
	{
		dist = dist + delta[i];
		ssf = this->man_lp_matrix_h.real_ss_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,dist);
		dist = dist - delta[i];

		dist = dist - delta[i];
		ssb = this->man_lp_matrix_h.real_ss_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,dist);
		dist = dist + delta[i];

		double fdm = (ssf-ssb)/2./delta[i];
		printf("Fdm ss - Z : %20.12lf\n",fdm);

		double anal = this->man_lp_matrix_h.real_ss_grad_z_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,dist);
		printf("Ana ss - Z : %20.12lf\n",anal);
		printf("dr  : %20.12lf\n",delta[i]);
		printf("Err : %20.12lf % \n",(anal-fdm)/anal*100.);
	}

	cout << endl;
	for(int i=0;i<6;i++)
	{
		dist = dist + delta[i];
		ssf = this->man_lp_matrix_h.real_sz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,dist);
		dist = dist - delta[i];

		dist = dist - delta[i];
		ssb = this->man_lp_matrix_h.real_sz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,dist);
		dist = dist + delta[i];

		double fdm = (ssf-ssb)/2./delta[i];
		printf("Fdm sz - Z : %20.12lf\n",fdm);

		double anal = this->man_lp_matrix_h.real_sz_grad_z_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,dist);
		printf("Ana sz - Z : %20.12lf\n",anal);
		printf("dr  : %20.12lf\n",delta[i]);
		printf("Err : %20.12lf % \n",(anal-fdm)/anal*100.);
	}

	cout << endl;
	for(int i=0;i<6;i++)
	{
		dist = dist + delta[i];
		ssf = this->man_lp_matrix_h.real_xx_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,dist);
		dist = dist - delta[i];

		dist = dist - delta[i];
		ssb = this->man_lp_matrix_h.real_xx_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,dist);
		dist = dist + delta[i];

		double fdm = (ssf-ssb)/2./delta[i];
		printf("Fdm xx - Z : %20.12lf\n",fdm);

		double anal = this->man_lp_matrix_h.real_xx_grad_z_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,dist);
		printf("Ana xx - Z : %20.12lf\n",anal);
		printf("dr  : %20.12lf\n",delta[i]);
		printf("Err : %20.12lf % \n",(anal-fdm)/anal*100.);
	}

	cout << endl;
	for(int i=0;i<6;i++)
	{
		dist = dist + delta[i];
		ssf = this->man_lp_matrix_h.real_zz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,dist);
		dist = dist - delta[i];

		dist = dist - delta[i];
		ssb = this->man_lp_matrix_h.real_zz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,dist);
		dist = dist + delta[i];

		double fdm = (ssf-ssb)/2./delta[i];
		printf("Fdm zz - Z : %20.12lf\n",fdm);

		double anal = this->man_lp_matrix_h.real_zz_grad_z_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,dist);
		printf("Ana zz - Z : %20.12lf\n",anal);
		printf("dr  : %20.12lf\n",delta[i]);
		printf("Err : %20.12lf % \n",(anal-fdm)/anal*100.);
	}
	exit(1);
#endif

	// Method Actual...
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

//	Optimise Periodic Summation Workload
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


////	////	////	////	////	////	////	////	////	////	////	////	////

////	Coulomb Interaction ( Periodic Summation )

////	////	////	////	////	////	////	////	////	////	////	////	////

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



////	////	////	////	////	////	////

////	LonePair_Member_Functions

////	////	////	////	////	////	////

void Manager::InitialiseLonePairEnergy( Cell& C )
{
	for(int i=0;i<C.NumberOfAtoms;i++)
	{
		if( C.AtomList[i]->type == "lone" )
		{	
			// Setting Temporal h_matrix zeroes
			static_cast<LonePair*>(C.AtomList[i])->lp_h_matrix_tmp.setZero();

			// Setting Onsite LonePair Model Parameter Lambda
			static_cast<LonePair*>(C.AtomList[i])->lp_h_matrix_tmp(1,1) = static_cast<LonePair*>(C.AtomList[i])->lp_lambda;
			static_cast<LonePair*>(C.AtomList[i])->lp_h_matrix_tmp(2,2) = static_cast<LonePair*>(C.AtomList[i])->lp_lambda;
			static_cast<LonePair*>(C.AtomList[i])->lp_h_matrix_tmp(3,3) = static_cast<LonePair*>(C.AtomList[i])->lp_lambda;
		}
	}
} 

void Manager::InitialiseLonePairDerivative( Cell& C ) {}
void Manager::InitialiseSCF() { this->man_vec.clear(); }	// Clear man_vec "Manager_Vector"

bool Manager::IsSCFDone( const double tol )			// Check If SCF Converged
{
	if( this->man_vec.size() < 2 )	// i.e., IF THIS IS THE 'FIRST SCF CYCLE'
	{	return false;
	}
	else	// IF IS THE CYCLES AFTHER THE FIRST
	{	if( fabs(this->man_vec[this->man_vec.size()-1] - this->man_vec[this->man_vec.size()-2]) > tol ) { return false; } // IF THE RECENT ENERGY PAIR DIFFERENCE IS LESS THAN THE TOLERANCE
		else
		{ 	//this->man_vec.claer();
			return true; 
		}
	}
}

void Manager::GetLonePairGroundState( Cell& C )	// Including Matrix Diagonalisaion + SetGroundState Index
{
	// This Function Assumes "LonePair::Eigen::Matrix4d lp_h_matrix_tmp" is Ready to be diagonalised

	std::vector<double> v(4);	// SPACE FOR HOLDING EIGEN VALUES
	double lp_scf_sum = 0.;		// TEMOPORALILY HOLDING GROUND STATE ENERGY SUM

	for(int i=0;i<C.NumberOfAtoms;i++)
	{
		if( C.AtomList[i]->type == "lone" )
		{	
			static_cast<LonePair*>(C.AtomList[i])->lp_h_matrix = static_cast<LonePair*>(C.AtomList[i])->lp_h_matrix_tmp;
			// Copy lp_h_matrix_tmp -> (into) lp_h_matrix ... dialgonalisation target
			static_cast<LonePair*>(C.AtomList[i])->lp_eigensolver.compute(static_cast<LonePair*>(C.AtomList[i])->lp_h_matrix,true);
			// Diagonalise LonePair H matrix, compute_evec='true'
		

			// Get LonePair GroundState Index
			v[0] = static_cast<LonePair*>(C.AtomList[i])->lp_eigensolver.eigenvalues()(0).real();
			v[1] = static_cast<LonePair*>(C.AtomList[i])->lp_eigensolver.eigenvalues()(1).real();
			v[2] = static_cast<LonePair*>(C.AtomList[i])->lp_eigensolver.eigenvalues()(2).real();
			v[3] = static_cast<LonePair*>(C.AtomList[i])->lp_eigensolver.eigenvalues()(3).real();
			
			static_cast<LonePair*>(C.AtomList[i])->lp_gs_index = std::min_element(v.begin(),v.end()) - v.begin();
			
			lp_scf_sum += static_cast<LonePair*>(C.AtomList[i])->lp_eigensolver.eigenvalues()(static_cast<LonePair*>(C.AtomList[i])->lp_gs_index).real();
		}
	}
	this->man_vec.push_back(lp_scf_sum);	// Logging CycSum
}

////	////	////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////

////	LonePair_Energy

void Manager::CoulombLonePairReal( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector, const bool is_first_scf )
{
	double Qi,Qj;
	Eigen::Vector3d Rij;
        // TransVector = h*a + k*b + l*c
        // Rij         = Ai.r - Aj.r - TransVector;

        if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {
		if( is_first_scf == true )
		{
			// LonePair Core - Core
			Qi  = C.AtomList[i]->charge;	// Get Lone CoreCharge
			Qj  = C.AtomList[j]->charge;	// Get CoreCharge
			Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
			C.mono_real_energy += 0.5*(Qi*Qj)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;
		}
		
		// LonePair - Core
		Qi  = static_cast<LonePair*>(C.AtomList[i])->lp_charge;
		Qj  = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;	// Ri(LonePairCore) - Rj(Core) - TransVector;
		// Calculate ... Matrix Elements ... add into : static_cast<LonePair*>(C.AtomList[i])->lp_h_matrix_tmp(a,b)
 

	}
	
	if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "lone" )
	{
		if( is_first_scf == true )
		{
			// Core - LonePair Core
			Qi  = C.AtomList[i]->charge;	// Get CoreCharge
			Qj  = C.AtomList[j]->charge;	// Get Lone CoreCharge
			Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
			C.mono_real_energy += 0.5*(Qi*Qj)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;
		}

	}

	if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "shel" )
	{
		if( is_first_scf == true )
		{
			// LonePair Core - Shell Core
			Qi  = C.AtomList[i]->charge;	// Get Lone CoreCharge
			Qj  = C.AtomList[j]->charge;	// Get Shel CoreCharge
			Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
			C.mono_real_energy += 0.5*(Qi*Qj)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;

			// LonePair Core - Shell Shell
			Qi  = C.AtomList[i]->charge;
			Qj  = static_cast<Shell*>(C.AtomList[j])->shel_charge;
			Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart - TransVector;
			C.mono_real_energy += 0.5*(Qi*Qj)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;
		}

	}

	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "lone" )
	{
		if( is_first_scf == true )
		{
			// Shell Core - LonePair Core
			Qi  = C.AtomList[i]->charge;
			Qj  = C.AtomList[j]->charge;
			Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
			C.mono_real_energy += 0.5*(Qi*Qj)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;

			// Shell Shell - LonePair Core
			Qi  = static_cast<Shell*>(C.AtomList[i])->shel_charge;
			Qj  = C.AtomList[j]->charge;
			Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - C.AtomList[j]->cart - TransVector;
			C.mono_real_energy += 0.5*(Qi*Qj)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;
		}

	}

	if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "lone" )
	{
		if( is_first_scf == true )
		{
			// LonePair Core - LonePair Core
			Qi  = C.AtomList[i]->charge;
			Qj  = C.AtomList[j]->charge;
			Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
			C.mono_real_energy += 0.5*(Qi*Qj)/Rij.norm() * erfc(Rij.norm()/C.sigma) * C.TO_EV;
		}
	}
}

void Manager::CoulombLonePairSelf( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector, const bool is_first_scf )
{
	double Qi,Qj;
	Eigen::Vector3d Rij;
        // TransVector = h*a + k*b + l*c
        // Rij         = Ai.r - Aj.r - TransVector;

	if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "lone" )	// SelfEnergy by LonePair Cores
	{
		if( is_first_scf == true )
		{
			Qi  = C.AtomList[i]->charge;
			Qj  = C.AtomList[j]->charge;
			C.mono_reci_self_energy += -0.5*(Qi*Qj)*2./C.sigma/sqrt(M_PI) * C.TO_EV;
		}
	}
}

void Manager::CoulombLonePairReci( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector, const bool is_first_scf )
{
	double Qi,Qj;
	double g_norm = TransVector.norm();
	double g_sqr  = g_norm*g_norm;
	Eigen::Vector3d Rij;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r;

        if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {
		if( is_first_scf == true )
		{
			// LonePair Core - Core
			Qi  = C.AtomList[i]->charge;	// Get Lone CoreCharge
			Qj  = C.AtomList[j]->charge;	// Get CoreCharge
			Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;
			C.mono_reci_energy += (2.*M_PI/C.volume)*(Qi*Qj)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;
		}
	}
	
	if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "lone" )
	{
		if( is_first_scf == true )
		{
			// Core - LonePair Core
			Qi  = C.AtomList[i]->charge;	// Get CoreCharge
			Qj  = C.AtomList[j]->charge;	// Get Lone CoreCharge
			Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;
			C.mono_reci_energy += (2.*M_PI/C.volume)*(Qi*Qj)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;
		}

	}

	if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "shel" )
	{
		if( is_first_scf == true )
		{
			// LonePair Core - Shell Core
			Qi  = C.AtomList[i]->charge;	// Get Lone CoreCharge
			Qj  = C.AtomList[j]->charge;	// Get Shel CoreCharge
			Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;
			C.mono_reci_energy += (2.*M_PI/C.volume)*(Qi*Qj)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;

			// LonePair Core - Shell Shell
			Qi  = C.AtomList[i]->charge;
			Qj  = static_cast<Shell*>(C.AtomList[j])->shel_charge;
			Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart;
			C.mono_reci_energy += (2.*M_PI/C.volume)*(Qi*Qj)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;
		}

	}

	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "lone" )
	{
		if( is_first_scf == true )
		{
			// Shell Core - LonePair Core
			Qi  = C.AtomList[i]->charge;
			Qj  = C.AtomList[j]->charge;
			Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;
			C.mono_reci_energy += (2.*M_PI/C.volume)*(Qi*Qj)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;

			// Shell Shell - LonePair Core
			Qi  = static_cast<Shell*>(C.AtomList[i])->shel_charge;
			Qj  = C.AtomList[j]->charge;
			Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - C.AtomList[j]->cart;
			C.mono_reci_energy += (2.*M_PI/C.volume)*(Qi*Qj)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;
		}

	}

	if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "lone" )
	{
		if( is_first_scf == true )
		{
			// LonePair Core - LonePair Core
			Qi  = C.AtomList[i]->charge;
			Qj  = C.AtomList[j]->charge;
			Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;
			C.mono_reci_energy += (2.*M_PI/C.volume)*(Qi*Qj)*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * cos( TransVector.adjoint()*Rij ) * C.TO_EV;
		}
	}
}

////	////	////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////

////	LonePair_Derivative

void Manager::CoulombLonePairDerivativeReal( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector )
{
	double Qi,Qj;
	double r_norm,r_sqr;
	Eigen::Vector3d Rij;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r - TransVector;
	double intact;

        if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {	
		// LonePair Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.AtomList[i]->cart_gd += intact*Rij;
		C.AtomList[j]->cart_gd -= intact*Rij;
        }       
        
	if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "lone" ) 
        {
		// Core - LonePair Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;
	
		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.AtomList[i]->cart_gd += intact*Rij;
		C.AtomList[j]->cart_gd -= intact*Rij;
        }       
        
	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "lone" ) 
        {
		// Shell Core - LonePair Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;
	
		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.AtomList[i]->cart_gd += intact*Rij;
		C.AtomList[j]->cart_gd -= intact*Rij;

		// Shell Shell - LonePair Core
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = C.AtomList[j]->charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		static_cast<Shell*>(C.AtomList[i])->shel_cart_gd += intact*Rij;
		C.AtomList[j]->cart_gd -= intact*Rij;
        }       
        
	if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "shel" ) 
        {
		// LonePair Core - Shell Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;
	
		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.AtomList[i]->cart_gd += intact*Rij;
		C.AtomList[j]->cart_gd -= intact*Rij;
		
		// LonePair Core - Shell Shell
		Qi = C.AtomList[i]->charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.AtomList[i]->cart_gd += intact*Rij;
		static_cast<Shell*>(C.AtomList[j])->shel_cart_gd -= intact*Rij;
        }

        if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "lone" )    // Handling Core - Core (i.e., charge charge interaction);
        {	
		// LonePair Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.AtomList[i]->cart_gd += intact*Rij;
		C.AtomList[j]->cart_gd -= intact*Rij;
        }       

}

void Manager::CoulombLonePairDerivativeSelf( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector )
{
	double Qi,Qj;
	double r_norm,r_sqr;
	Eigen::Vector3d Rij;
	
	double intact;

	if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "lone" ) 
        {

	}
}       

void Manager::CoulombLonePairDerivativeReci( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector )
{
	double Qi,Qj;
	double g_norm = TransVector.norm();
	double g_sqr  = g_norm*g_norm;
	Eigen::Vector3d Rij;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r;
	double intact;
        
        if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {	
		// LonePair Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;

		intact = C.TO_EV*((2.*M_PI)/C.volume)*Qi*Qj*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);

		C.AtomList[i]->cart_gd += intact * TransVector;
		C.AtomList[j]->cart_gd -= intact * TransVector;
        }       
        
	if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "lone" ) 
        {	
		// Core - LonePair Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;

		intact = C.TO_EV*((2.*M_PI)/C.volume)*Qi*Qj*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);

		C.AtomList[i]->cart_gd += intact * TransVector;
		C.AtomList[j]->cart_gd -= intact * TransVector;
        }       
        
	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "lone" ) 
        {
		// Shell Core - LonePair Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;

		intact = C.TO_EV*((2.*M_PI)/C.volume)*Qi*Qj*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);

		C.AtomList[i]->cart_gd += intact * TransVector;
		C.AtomList[j]->cart_gd -= intact * TransVector;

		// Shell Shell - LonePair Core
		Qi = static_cast<Shell*>(C.AtomList[i])->shel_charge;
		Qj = C.AtomList[j]->charge;
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - C.AtomList[j]->cart;
	
		intact = C.TO_EV*((2.*M_PI)/C.volume)*Qi*Qj*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);

		static_cast<Shell*>(C.AtomList[i])->shel_cart_gd += intact * TransVector;
		C.AtomList[j]->cart_gd -= intact * TransVector;
        }       
        
	if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "shel" ) 
        {
		// LonePair Core - Shell Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;

		intact = C.TO_EV*((2.*M_PI)/C.volume)*Qi*Qj*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);

		C.AtomList[i]->cart_gd += intact * TransVector;
		C.AtomList[j]->cart_gd -= intact * TransVector;

		// LonePair Core - Shell Shell
		Qi = C.AtomList[i]->charge;
		Qj = static_cast<Shell*>(C.AtomList[j])->shel_charge;
		Rij = C.AtomList[i]->cart - static_cast<Shell*>(C.AtomList[j])->shel_cart;

		intact = C.TO_EV*((2.*M_PI)/C.volume)*Qi*Qj*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);

		C.AtomList[i]->cart_gd += intact * TransVector;
		static_cast<Shell*>(C.AtomList[j])->shel_cart_gd -= intact * TransVector;
        }       

        if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "lone" )    // Handling Core - Core (i.e., charge charge interaction);
        {	
		// LonePair Core - Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart;

		intact = C.TO_EV*((2.*M_PI)/C.volume)*Qi*Qj*exp(-0.25*C.sigma*C.sigma*g_sqr)/g_sqr * -sin(TransVector.adjoint()*Rij);

		C.AtomList[i]->cart_gd += intact * TransVector;
		C.AtomList[j]->cart_gd -= intact * TransVector;
        }       
}       

////	////	////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////    ////

////	LonePair_StrainDerivative

void Manager::StrainLonePairDerivativeReal( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector )
{
	double Qi,Qj;
	double r_norm,r_sqr;
	Eigen::Vector3d Rij;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r - TransVector;
	double intact;
        
        if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {
		// LonePair Core - Core
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
        
	if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "shel" ) 
        {	
		// LonePair Core - Shell Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.lattice_sd(0,0) += intact * Rij(0) * Rij(0);	C.lattice_sd(0,1) += intact * Rij(0) * Rij(1);	C.lattice_sd(0,2) += intact * Rij(0) * Rij(2);
								C.lattice_sd(1,1) += intact * Rij(1) * Rij(1);	C.lattice_sd(1,2) += intact * Rij(1) * Rij(2);
														C.lattice_sd(2,2) += intact * Rij(2) * Rij(2);
		// LonePair Core - Shell Shell
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
        
	if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "lone" ) 
        {	
		// Core - LonePair Core
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
        
	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "lone" ) 
        {
		// Shell Core - LonePair Core
		Qi = C.AtomList[i]->charge;
		Qj = C.AtomList[j]->charge;
		Rij = C.AtomList[i]->cart - C.AtomList[j]->cart - TransVector;
		r_norm = Rij.norm();
		r_sqr  = r_norm*r_norm;

		intact = C.TO_EV*(0.5*Qi*Qj)*((-2./C.sigma/sqrt(M_PI))*(exp(-r_sqr/C.sigma/C.sigma)/r_norm)-(erfc(r_norm/C.sigma)/r_sqr))/r_norm;

		C.lattice_sd(0,0) += intact * Rij(0) * Rij(0);	C.lattice_sd(0,1) += intact * Rij(0) * Rij(1);	C.lattice_sd(0,2) += intact * Rij(0) * Rij(2);
								C.lattice_sd(1,1) += intact * Rij(1) * Rij(1);	C.lattice_sd(1,2) += intact * Rij(1) * Rij(2);
														C.lattice_sd(2,2) += intact * Rij(2) * Rij(2);
		// Shell Shell - LonePair Core
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

        if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "lone" )    // Handling Core - Core (i.e., charge charge interaction);
        {
		// LonePair Core - Core
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
}       

void Manager::StrainLonePairDerivativeSelf( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector )
{
	double Qi,Qj;
	double r_norm,r_sqr;
	Eigen::Vector3d Rij;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r - TransVector;
	double intact;

        if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "lone" )    // Handling Core - Core (i.e., charge charge interaction);
	{

	}
}

void Manager::StrainLonePairDerivativeReci( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector )
{
	double Qi,Qj;
	double g_norm = TransVector.norm();
	double g_sqr  = g_norm*g_norm;
	Eigen::Vector3d Rij;
        // TransVector(G) = 2pi h*u + 2pi k*v + 2pi l*w;
        // Rij            = Ai.r - Aj.r;
	double intact[4];

        if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "core" )    // Handling Core - Core (i.e., charge charge interaction);
        {
		// LonePair Core - Core
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
        
	if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "shel" ) 
        {
		// LonePair Core - Shell Core
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

		// LonePair Core - Shell Shell
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
        
	if( C.AtomList[i]->type == "core" && C.AtomList[j]->type == "lone" ) 
        {
		// Core - LonePair Core
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
        
	if( C.AtomList[i]->type == "shel" && C.AtomList[j]->type == "lone" ) 
        {
		// Shell Core - LonePair Core
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

		// Shell Shell - LonePair Core
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

        if( C.AtomList[i]->type == "lone" && C.AtomList[j]->type == "lone" )    // Handling Core - Core (i.e., charge charge interaction);
        {
		// LonePair Core - Core
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
}
