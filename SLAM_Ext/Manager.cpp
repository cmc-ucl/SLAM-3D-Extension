#include "Manager.hpp"

//#define PRINT_REAL_SPACE_D2

#define DEV_G_SPACE

// Implement RealSpace Integrators - Input ... LonePair* / Vector to a species / sigma / IndexLonePair* / IndexSpecies

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
}

void Manager::set_h_matrix_real_pc_derivative( /* IN/RES OUT */ LonePair* lp, const Eigen::Vector3d& R, const double sig, const int lp_i, const int pc_i )
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
	
	// 1. Compute first derivative integrals in a local symmetry
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

	// 2. Using the local elements; compute equivalent elements (in the global) in the local reference frame
	// note : h_tmp_*_ws are in the local reference frame, their x'/y'/z' element (local) has to be inversed to x/y/z (global) 
	h_tmp_x_ws = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp_x * this->man_lp_matrix_h.transform_matrix;
	h_tmp_y_ws = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp_y * this->man_lp_matrix_h.transform_matrix;
	h_tmp_z_ws = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp_z * this->man_lp_matrix_h.transform_matrix;

	// 3. Transform back to the global reference frame
	for(int i=0;i<4;i++)
	{	for(int j=0;j<4;j++)
		{	v_loc << h_tmp_x_ws(i,j), h_tmp_y_ws(i,j), h_tmp_z_ws(i,j);
			v_glo = this->man_lp_matrix_h.transform_matrix_shorthand.transpose() * v_loc;
			this->real_lp_h_pc_x[lp_i][pc_i](i,j) =  v_glo(0);
			this->real_lp_h_pc_y[lp_i][pc_i](i,j) =  v_glo(1);
			this->real_lp_h_pc_z[lp_i][pc_i](i,j) =  v_glo(2);
		}
	}
}

void Manager::set_h_matrix_real_pc_derivative2( LonePair* lp, const Eigen::Vector3d& R, const double sig, const int lp_i, const int pc_i )
{
	this->man_lp_matrix_h.GetTransformationMatrix(R);	// get Transformation matrix ... saved : Eigen::Matrix4d this->man_lp_matrix_h.transform_matrix;
	Eigen::Matrix4d h_tmp_d2[3][3];
	Eigen::Matrix4d h_d2_ws[3][3];
	Eigen::Matrix3d m_loc,m_glo;

	for(int i=0;i<3;i++) { for(int j=0;j<3;j++) { h_tmp_d2[i][j].setZero(); h_d2_ws[i][j].setZero(); }}	// Initialise workspace

	// [0][0] xx h_matrix loc
	h_tmp_d2[0][0](0,0) = this->man_lp_matrix_h.real_ss_grad2_xx_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());	// [0][0] - XX, (0,0) ss
	h_tmp_d2[0][0](0,3) = this->man_lp_matrix_h.real_sz_grad2_xx_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());	// [0][0] - XX, (0,3) sz
	h_tmp_d2[0][0](3,0) = h_tmp_d2[0][0](0,3);
	h_tmp_d2[0][0](1,1) = this->man_lp_matrix_h.real_xx_grad2_xx_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());	// [0][0] - XX, (1,1) xx
	h_tmp_d2[0][0](2,2) = this->man_lp_matrix_h.real_yy_grad2_xx_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());	// [0][0] - XX, (2,2) yy
	h_tmp_d2[0][0](3,3) = this->man_lp_matrix_h.real_zz_grad2_xx_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());	// [0][0] - XX, (3,3) zz

	// [0][1] xy h_matrix loc
	h_tmp_d2[0][1](1,2) = this->man_lp_matrix_h.real_xy_grad2_xy_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());	// [0][1] - XY, (1,2) xy
	h_tmp_d2[0][1](2,1) = h_tmp_d2[0][1](1,2);

	// [0][2] xz h_matrix loc
	h_tmp_d2[0][2](0,1) = this->man_lp_matrix_h.real_sx_grad2_xz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());	// [0][2] - XZ, (0,1) sx
	h_tmp_d2[0][2](1,0) = h_tmp_d2[0][2](0,1);
	h_tmp_d2[0][2](1,3) = this->man_lp_matrix_h.real_xz_grad2_xz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());	// [0][2] - XZ, (1,3) xz
	h_tmp_d2[0][2](3,1) = h_tmp_d2[0][2](1,3);

	// [2][2] zz h_matix loc
	h_tmp_d2[2][2](0,0) = this->man_lp_matrix_h.real_ss_grad2_zz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());	// [2][2] - ZZ, (0,0) ss
	h_tmp_d2[2][2](0,3) = this->man_lp_matrix_h.real_sz_grad2_zz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());	// [2][2] - ZZ, (0,3) sz
	h_tmp_d2[2][2](3,0) = h_tmp_d2[2][2](0,3);
	h_tmp_d2[2][2](1,1) = this->man_lp_matrix_h.real_xx_grad2_zz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());	// [2][2] - ZZ, (0,3) xx
	h_tmp_d2[2][2](2,2) = h_tmp_d2[2][2](1,1);
	h_tmp_d2[2][2](3,3) = this->man_lp_matrix_h.real_zz_grad2_zz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());	// [2][2] - ZZ, (3,3) zz

	// Below here, done by the orbital symmetry
	// [1][0] yx - [0][1] xy
	h_tmp_d2[1][0] = h_tmp_d2[0][1];
	// [1][1] yy ~ [0][0] xx
	h_tmp_d2[1][1] = h_tmp_d2[0][0];
	h_tmp_d2[1][1](1,1) = h_tmp_d2[0][0](2,2);	// YY - xx = XX - yy
	h_tmp_d2[1][1](2,2) = h_tmp_d2[0][0](1,1);	// YY - yy = XX - xx
	// [1][2] yz ~ [0][2] xz
	h_tmp_d2[1][2](0,2) = h_tmp_d2[1][2](2,0) = h_tmp_d2[0][2](0,1);	// YZ - sy = XZ - sx
	h_tmp_d2[1][2](2,3) = h_tmp_d2[1][2](3,2) = h_tmp_d2[0][2](1,3);	// YZ - yz = XZ - xz
	// [2][0] zx - [0][2] xz
	h_tmp_d2[2][0] = h_tmp_d2[0][2];
	// [2][1] zy - [1][2] yz
	h_tmp_d2[2][1] = h_tmp_d2[1][2];
	//// End of Local element Set

	// 2. Using the local elements; compute equivalent elements (in the global) in the local reference frame
	for(int i=0;i<3;i++){ for(int j=0;j<3;j++){ h_d2_ws[i][j] = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp_d2[i][j] * this->man_lp_matrix_h.transform_matrix; }}

	// 3. Transform back to the global reference frame - i,j refer to basis functions
	for(int i=0;i<4;i++)
	{	for(int j=0;j<4;j++)
		{	
			m_loc << h_d2_ws[0][0](i,j), h_d2_ws[0][1](i,j), h_d2_ws[0][2](i,j),
				 h_d2_ws[1][0](i,j), h_d2_ws[1][1](i,j), h_d2_ws[1][2](i,j),
				 h_d2_ws[2][0](i,j), h_d2_ws[2][1](i,j), h_d2_ws[2][2](i,j);

			m_glo = this->man_lp_matrix_h.transform_matrix_shorthand.transpose() * m_loc * this->man_lp_matrix_h.transform_matrix_shorthand;

			this->real_lp_h_lp_xx[lp_i][pc_i](i,j) = m_glo(0,0); this->real_lp_h_lp_xy[lp_i][pc_i](i,j) = m_glo(0,1); this->real_lp_h_lp_xz[lp_i][pc_i](i,j) = m_glo(0,2);
			this->real_lp_h_lp_yx[lp_i][pc_i](i,j) = m_glo(1,0); this->real_lp_h_lp_yy[lp_i][pc_i](i,j) = m_glo(1,1); this->real_lp_h_lp_yz[lp_i][pc_i](i,j) = m_glo(1,2);
			this->real_lp_h_lp_zx[lp_i][pc_i](i,j) = m_glo(2,0); this->real_lp_h_lp_zy[lp_i][pc_i](i,j) = m_glo(2,1); this->real_lp_h_lp_zz[lp_i][pc_i](i,j) = m_glo(2,2);
		}
	}
}

const Eigen::Matrix4d& Manager::set_h_matrix_reci_cos( /* IN/RES OUT */ LonePair* lp, const Eigen::Vector3d& G, const int lp_i, const int pc_i )
{
	this->man_lp_matrix_h.GetTransformationMatrix(G);	// get Transformation matrix ... saved : Eigen::Matrix4d this->man_lp_matrix_h.transform_matrix;
	Eigen::Matrix4d h_tmp;					// calculating local h_matrix WS
	h_tmp.setZero();					// initialise

	// evalulation block
	h_tmp(0,0) = this->man_lp_matrix_h.reci_ss_cos(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,G.norm());
	h_tmp(1,1) = this->man_lp_matrix_h.reci_xx_cos(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,G.norm());
	h_tmp(2,2) = h_tmp(1,1);
	h_tmp(3,3) = this->man_lp_matrix_h.reci_zz_cos(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,G.norm());

	this->reci_lp_h_pc[lp_i][pc_i] = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp * this->man_lp_matrix_h.transform_matrix;	// inverse transformation

	return this->reci_lp_h_pc[lp_i][pc_i];
}


void Manager::set_h_matrix_reci_derivative_cos( /* IN/RES OUT */ LonePair* lp, const Eigen::Vector3d& G, const int lp_i, const int pc_i )
{
	this->man_lp_matrix_h.GetTransformationMatrix(G);	// get Transformation matrix ... saved : Eigen::Matrix4d this->man_lp_matrix_h.transform_matrix;
	Eigen::Matrix4d h_tmp_x, h_tmp_y, h_tmp_z;		// calculating local h_matrix WS
	Eigen::Matrix4d h_tmp_x_ws,h_tmp_y_ws,h_tmp_z_ws;
	Eigen::Vector3d v_loc, v_glo;

	h_tmp_x.setZero();					// initialise
	h_tmp_y.setZero();					// initialise
	h_tmp_z.setZero();					// initialise
	h_tmp_x_ws.setZero();					// initialise
	h_tmp_y_ws.setZero();					// initialise
	h_tmp_z_ws.setZero();					// initialise
	
	// 1. Compute first derivative integrals in a local symmetry
	h_tmp_x(1,3) = this->man_lp_matrix_h.reci_xz_grad_gx_cos(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,G.norm());
	h_tmp_x(3,1) = h_tmp_x(1,3);

	h_tmp_y(2,3) = h_tmp_y(3,2) = h_tmp_x(1,3);	// y-yz = x-sz
	
	h_tmp_z(0,0) = this->man_lp_matrix_h.reci_ss_grad_gz_cos(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,G.norm());
	h_tmp_z(1,1) = this->man_lp_matrix_h.reci_xx_grad_gz_cos(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,G.norm());
	h_tmp_z(2,2) = h_tmp_z(1,1);
	h_tmp_z(3,3) = this->man_lp_matrix_h.reci_zz_grad_gz_cos(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,G.norm());

	// 2. Using the local elements; compute equivalent elements (in the global) in the local reference frame
	// note : h_tmp_*_ws are in the local reference frame, their x'/y'/z' element (local) has to be inversed to x/y/z (global) 
	h_tmp_x_ws = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp_x * this->man_lp_matrix_h.transform_matrix;
	h_tmp_y_ws = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp_y * this->man_lp_matrix_h.transform_matrix;
	h_tmp_z_ws = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp_z * this->man_lp_matrix_h.transform_matrix;

	// 3. Transform back to the global reference frame
	for(int i=0;i<4;i++)
	{	for(int j=0;j<4;j++)
		{	v_loc << h_tmp_x_ws(i,j), h_tmp_y_ws(i,j), h_tmp_z_ws(i,j);
			v_glo = this->man_lp_matrix_h.transform_matrix_shorthand.transpose() * v_loc;
			this->reci_lp_h_pc_x[lp_i][pc_i](i,j) =  v_glo(0);
			this->reci_lp_h_pc_y[lp_i][pc_i](i,j) =  v_glo(1);
			this->reci_lp_h_pc_z[lp_i][pc_i](i,j) =  v_glo(2);
		}
	}
}


const Eigen::Matrix4d& Manager::set_h_matrix_reci_sin( /* IN/RES OUT */ LonePair* lp, const Eigen::Vector3d& G, const int lp_i, const int pc_i )
{
	this->man_lp_matrix_h.GetTransformationMatrix(G);	// get Transformation matrix ... saved : Eigen::Matrix4d this->man_lp_matrix_h.transform_matrix;
	Eigen::Matrix4d h_tmp;					// calculating local h_matrix WS
	h_tmp.setZero();					// initialise

	// evalulation block
	h_tmp(0,3) = this->man_lp_matrix_h.reci_sz_sin(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,G.norm());
	h_tmp(3,0) = h_tmp(0,3);

	this->reci_lp_h_pc[lp_i][pc_i] = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp * this->man_lp_matrix_h.transform_matrix;	// inverse transformation

	return this->reci_lp_h_pc[lp_i][pc_i];
}


void Manager::set_h_matrix_reci_derivative_sin( /* IN/RES OUT */ LonePair* lp, const Eigen::Vector3d& G, const int lp_i, const int pc_i )
{
	this->man_lp_matrix_h.GetTransformationMatrix(G);	// get Transformation matrix ... saved : Eigen::Matrix4d this->man_lp_matrix_h.transform_matrix;
	Eigen::Matrix4d h_tmp_x, h_tmp_y, h_tmp_z;		// calculating local h_matrix WS
	Eigen::Matrix4d h_tmp_x_ws,h_tmp_y_ws,h_tmp_z_ws;
	Eigen::Vector3d v_loc, v_glo;

	h_tmp_x.setZero();					// initialise
	h_tmp_y.setZero();					// initialise
	h_tmp_z.setZero();					// initialise
	h_tmp_x_ws.setZero();					// initialise
	h_tmp_y_ws.setZero();					// initialise
	h_tmp_z_ws.setZero();					// initialise
	
	// 1. Compute first derivative integrals in a local symmetry
	h_tmp_x(0,1) = this->man_lp_matrix_h.reci_sx_grad_gx_sin(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,G.norm());
	h_tmp_x(1,0) = h_tmp_x(0,1);

	h_tmp_y(0,2) = h_tmp_y(2,0) = h_tmp_x(0,1);	// y-sy = x-sx
	
	h_tmp_z(0,3) = this->man_lp_matrix_h.reci_sz_grad_gz_sin(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,G.norm());
	h_tmp_z(3,0) = h_tmp_z(0,3);

	// 2. Using the local elements; compute equivalent elements (in the global) in the local reference frame
	// note : h_tmp_*_ws are in the local reference frame, their x'/y'/z' element (local) has to be inversed to x/y/z (global) 
	h_tmp_x_ws = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp_x * this->man_lp_matrix_h.transform_matrix;
	h_tmp_y_ws = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp_y * this->man_lp_matrix_h.transform_matrix;
	h_tmp_z_ws = this->man_lp_matrix_h.transform_matrix.transpose() * h_tmp_z * this->man_lp_matrix_h.transform_matrix;

	// 3. Transform back to the global reference frame
	for(int i=0;i<4;i++)
	{	for(int j=0;j<4;j++)
		{	v_loc << h_tmp_x_ws(i,j), h_tmp_y_ws(i,j), h_tmp_z_ws(i,j);
			v_glo = this->man_lp_matrix_h.transform_matrix_shorthand.transpose() * v_loc;
			this->reci_lp_h_pc_x[lp_i][pc_i](i,j) =  v_glo(0);
			this->reci_lp_h_pc_y[lp_i][pc_i](i,j) =  v_glo(1);
			this->reci_lp_h_pc_z[lp_i][pc_i](i,j) =  v_glo(2);
		}
	}
}


void Manager::InitialiseEnergy( Cell& C )
{


#ifdef DEV_G_SPACE
using std::cout, std::endl;
LonePair* lp = nullptr;
int lp_id;

// FDM CHECK VARS
double delta = 0.005;
double sig   = 1.85;
Eigen::Vector3d v;

for(int i=0;i<C.NumberOfAtoms;i++)
{	//cout << "index : " << i+1 << " / type : " << C.AtomList[i]->type << endl;
	if( C.AtomList[i]->type == "lone" )
	{	lp_id = i;
		lp    = static_cast<LonePair*>(C.AtomList[i]);
		break;
	}
}
// set a vector
v << 1,2,3;
//auto begin = std::chrono::system_clock::now();
//set_h_matrix_real_pc(lp,v,sig,lp_id,0);
//set_h_matrix_real_pc_derivative(lp,v,sig,lp_id,0);
//set_h_matrix_real_pc_derivative2(lp,v,sig,lp_id,0);
//auto end   = std::chrono::system_clock::now() - begin;
//auto time_s= std::chrono::duration<double>(end).count();
//cout << "Elapsed time (s) : lp - lp interaction" << endl;
//cout << time_s << " (s)" << endl;

double g = 0.05;

cout << "G Vector" << endl;
//h_tmp_d2[0][2](0,1) = this->man_lp_matrix_h.real_sx_grad2_xz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig,R.norm());	// [0][2] - XZ, (0,1) sx

double ss,xx,zz;
double x_xz;
double z_ss,z_xx,z_zz;
/*
for(int i=0;i<1000;i++)
{
	g = g*1.0123;

	ss = this->man_lp_matrix_h.reci_ss_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,g);
	xx = this->man_lp_matrix_h.reci_xx_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,g);
	zz = this->man_lp_matrix_h.reci_zz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,g);
	x_xz = this->man_lp_matrix_h.reci_xz_grad_gx_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,g);
	z_ss = this->man_lp_matrix_h.reci_ss_grad_gz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,g);
	z_xx = this->man_lp_matrix_h.reci_xx_grad_gz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,g);
	z_zz = this->man_lp_matrix_h.reci_zz_grad_gz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,g);


	printf("%12.8e\t%20.12e\t%20.12e\t%20.12e\t%20.12e\t%20.12e\t%20.12e\t%20.12e\n",
		g,ss,xx,zz,z_ss,z_xx,z_zz,x_xz);
	
	if( g > 15. )
		break;
}
*/


cout << "G FDM CHECK" << endl;
double ss_f,xx_f,zz_f;
double ss_b,xx_b,zz_b;
g = 25.124;
delta = 0.0001;
cout << "FDM Z Forward" << endl;
g = g + delta;
ss_f = this->man_lp_matrix_h.reci_ss_cos(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,g);
xx_f = this->man_lp_matrix_h.reci_xx_cos(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,g);
zz_f = this->man_lp_matrix_h.reci_zz_cos(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,g);
printf("%20.12e\t%20.12e\t%20.12e\n",ss_f,xx_f,zz_f);
g = g - delta;
cout << "FDM Z Backward" << endl;
g = g - delta;
ss_b = this->man_lp_matrix_h.reci_ss_cos(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,g);
xx_b = this->man_lp_matrix_h.reci_xx_cos(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,g);
zz_b = this->man_lp_matrix_h.reci_zz_cos(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,g);
printf("%20.12e\t%20.12e\t%20.12e\n",ss_b,xx_b,zz_b);
g = g + delta;
cout << "FDM Onsite" << endl;
z_ss = this->man_lp_matrix_h.reci_ss_grad_gz_cos(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,g);
z_xx = this->man_lp_matrix_h.reci_xx_grad_gz_cos(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,g);
z_zz = this->man_lp_matrix_h.reci_zz_grad_gz_cos(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,g);
printf("%20.12e\t%20.12e\t%20.12e\n",z_ss,z_xx,z_zz);
printf("%20.12e\t%20.12e\t%20.12e\n",(ss_f-ss_b)/2./delta,(xx_f-xx_b)/2./delta,(zz_f-zz_b)/2./delta);

cout << endl;
cout << "Transformation Check" << endl;
Eigen::Vector3d G;
G << -1.14,0.44,0.3;
cout << endl;
cout << "H Onsite" << endl;
this->set_h_matrix_reci_cos(lp,G,lp_id,0);
for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ printf("%20.12e\t",this->reci_lp_h_pc[lp_id][0](i,j)); } cout << endl;}
cout << "dH Onsite" << endl;
this->set_h_matrix_reci_derivative_cos(lp,G,lp_id,0);
cout << "d/dgx" << endl;
for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ printf("%20.12e\t",this->reci_lp_h_pc_x[lp_id][0](i,j)); } cout << endl;}
cout << "d/dgy" << endl;
for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ printf("%20.12e\t",this->reci_lp_h_pc_y[lp_id][0](i,j)); } cout << endl;}
cout << "d/dgz" << endl;
for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ printf("%20.12e\t",this->reci_lp_h_pc_z[lp_id][0](i,j)); } cout << endl;}

cout << "Initiate FDM G" << endl;
Eigen::Matrix4d gx,gy,gz;
gx.setZero();
gy.setZero();
gz.setZero();

// X
cout << "FDM gx" << endl;
G(0) = G(0) + delta;
gx = this->set_h_matrix_reci_cos(lp,G,lp_id,0);
G(0) = G(0) - delta;
G(0) = G(0) - delta;
gx = gx - this->set_h_matrix_reci_cos(lp,G,lp_id,0);
G(0) = G(0) + delta;
gx = gx/2./delta;
for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ printf("%20.12e\t",gx(i,j)); } cout << endl;}
// Y
cout << "FDM gy" << endl;
G(1) = G(1) + delta;
gy = this->set_h_matrix_reci_cos(lp,G,lp_id,0);
G(1) = G(1) - delta;
G(1) = G(1) - delta;
gy = gy - this->set_h_matrix_reci_cos(lp,G,lp_id,0);
G(1) = G(1) + delta;
gy = gy/2./delta;
for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ printf("%20.12e\t",gy(i,j)); } cout << endl;}
// Z
cout << "FDM gz" << endl;
G(2) = G(2) + delta;
gz = this->set_h_matrix_reci_cos(lp,G,lp_id,0);
G(2) = G(2) - delta;
G(2) = G(2) - delta;
gz = gz - this->set_h_matrix_reci_cos(lp,G,lp_id,0);
G(2) = G(2) + delta;
gz = gz/2./delta;
for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ printf("%20.12e\t",gz(i,j)); } cout << endl;}



cout << endl;
cout << "PosIntegral test" << endl;
double posint;
posint = this->man_lp_matrix_h.real_position_integral(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function);
printf("PosIntegral : %20.12lf\n",posint);

double lp_self_ss;
double lp_self_xx;
double lp_self_grad;

sig = 1.3214;

lp_self_ss = this->man_lp_matrix_h.reci_self_integral_ss(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig);
lp_self_xx = this->man_lp_matrix_h.reci_self_integral_xx(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig);
lp_self_grad = this->man_lp_matrix_h.reci_self_integral_sx_grad_x(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sig);

cout << "lp_self ss / xx / lp_self_grad " << endl;
printf("%18.12lf\t%18.12lf\t%18.12lf\n",lp_self_ss,lp_self_xx,lp_self_grad);

cout << "### Terminating G-Space Recipe Dev" << endl;

cout << endl;
cout << endl;
cout << endl;
cout << "Testing SinPart" << endl;

cout << "G FDM CHECK" << endl;
delta = 0.0001;
G << -3.14,4.44,-0.3;
sig = 2.4214;
cout << "G Vector on Test: ";
printf("%18.12lf\t%18.12lf\t%18.12lf\n",G(0),G(1),G(2));
cout << "GNorm :           " << G.norm() << endl;

cout << endl;
cout << "H Onsite" << endl;
this->set_h_matrix_reci_sin(lp,G,lp_id,0);
for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ printf("%20.12e\t",this->reci_lp_h_pc[lp_id][0](i,j)); } cout << endl;}
cout << "dH Onsite" << endl;
this->set_h_matrix_reci_derivative_sin(lp,G,lp_id,0);
cout << "d/dgx" << endl;
for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ printf("%20.12e\t",this->reci_lp_h_pc_x[lp_id][0](i,j)); } cout << endl;}
cout << "d/dgy" << endl;
for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ printf("%20.12e\t",this->reci_lp_h_pc_y[lp_id][0](i,j)); } cout << endl;}
cout << "d/dgz" << endl;
for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ printf("%20.12e\t",this->reci_lp_h_pc_z[lp_id][0](i,j)); } cout << endl;}

cout << "Initiate FDM G" << endl;
//Eigen::Matrix4d gx,gy,gz;
gx.setZero();
gy.setZero();
gz.setZero();

// X
cout << "FDM gx" << endl;
G(0) = G(0) + delta;
gx = this->set_h_matrix_reci_sin(lp,G,lp_id,0);
G(0) = G(0) - delta;
G(0) = G(0) - delta;
gx = gx - this->set_h_matrix_reci_sin(lp,G,lp_id,0);
G(0) = G(0) + delta;
gx = gx/2./delta;
for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ printf("%20.12e\t",gx(i,j)); } cout << endl;}
// Y
cout << "FDM gy" << endl;
G(1) = G(1) + delta;
gy = this->set_h_matrix_reci_sin(lp,G,lp_id,0);
G(1) = G(1) - delta;
G(1) = G(1) - delta;
gy = gy - this->set_h_matrix_reci_sin(lp,G,lp_id,0);
G(1) = G(1) + delta;
gy = gy/2./delta;
for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ printf("%20.12e\t",gy(i,j)); } cout << endl;}
// Z
cout << "FDM gz" << endl;
G(2) = G(2) + delta;
gz = this->set_h_matrix_reci_sin(lp,G,lp_id,0);
G(2) = G(2) - delta;
G(2) = G(2) - delta;
gz = gz - this->set_h_matrix_reci_sin(lp,G,lp_id,0);
G(2) = G(2) + delta;
gz = gz/2./delta;
for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ printf("%20.12e\t",gz(i,j)); } cout << endl;}

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

/* Supportive Functions */
void Manager::support_h_matrix_real( const LonePair* lp, const double& sigma, const Eigen::Vector3d& Rij, /* workspace */ Eigen::Matrix4d& h_mat_ws, /* out */ Eigen::Matrix4d& h_mat_out )
{
	this->man_lp_matrix_h.GetTransformationMatrix(Rij);
	h_mat_ws.setZero();

	// Evaluation
	h_mat_ws(0,0) = this->man_lp_matrix_h.real_ss_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sigma,Rij.norm());
	h_mat_ws(0,3) = this->man_lp_matrix_h.real_sz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sigma,Rij.norm());
	h_mat_ws(3,0) = h_mat_ws(0,3);
	h_mat_ws(1,1) = this->man_lp_matrix_h.real_xx_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sigma,Rij.norm());
	h_mat_ws(2,2) = h_mat_ws(1,1);
	h_mat_ws(3,3) = this->man_lp_matrix_h.real_zz_pc(lp->lp_r,lp->lp_r_s_function,lp->lp_r_p_function,sigma,Rij.norm());

	// Inverse Transformation
	h_mat_out = this->man_lp_matrix_h.transform_matrix.transpose() * h_mat_ws * this->man_lp_matrix_h.transform_matrix;
}

void Manager::set_h_matrix_real( Cell& C, const int i, const int j, const Eigen::Vector3d& TransVector )
{
	const std::string type_i = C.AtomList[i]->type;
	const std::string type_j = C.AtomList[j]->type;
	double lp_cf[4];
	double factor;		// Multiplication Factor ... controls species charges + etcs.
	double partial_e;	// Partial energy when 'i' is non-LP species and 'j' is LP density
	Eigen::Vector3d Rij;

	/* if 'lone' comes to the place 'i'    : compute h matrix 

	   else (i.e., comes to the place 'j') : compute the interaction energy -> based on the previous charge density shape (i.e., LP lone pair eigenvectors)
	*/

	if( type_i == "lone" && type_j == "core" )	// LonePairD - Core 
	{
		LonePair* lp = static_cast<LonePair*>(C.AtomList[i]);

		// W.R.T Core
		Rij = ( C.AtomList[j]->cart + TransVector ) - C.AtomList[i]->cart;	// (Rj+T) - Ri ... Not using 'Ri - Rj - T' to get the right transformation // 'i' LPcore - 'j' core
		// Evaluation
		Manager::support_h_matrix_real( lp, C.sigma, Rij, this->man_matrix4d_ws[0], this->man_matrix4d_ws[1] );

		factor = 0.5 * lp->lp_charge * C.AtomList[j]->charge;			// Inverse Transformation .. mul 1/2 ... factor-out double counting
		this->LPC_H_Real[i][j][0] += factor * this->man_matrix4d_ws[1];
		/*
			H matrix element of 'i'th LP ? - should it be additive?
		*/
	}

	if( type_i == "core" && type_j == "lone" )	// LonePairD - Core
	{
		// Get 'j' lone pair cation & its eigenvectors of the ground-state
		LonePair* lp = static_cast<LonePair*>(C.AtomList[j]);
		lp_cf[0] = lp->lp_eigensolver.eigenvectors()(0,lp->lp_gs_index).real();	// s
		lp_cf[1] = lp->lp_eigensolver.eigenvectors()(1,lp->lp_gs_index).real();	// px
		lp_cf[2] = lp->lp_eigensolver.eigenvectors()(2,lp->lp_gs_index).real();	// py
		lp_cf[3] = lp->lp_eigensolver.eigenvectors()(3,lp->lp_gs_index).real();	// pz

		// W.R.T LP Density
		Rij = C.AtomList[i]->cart - ( C.AtomList[j]->cart + TransVector );	// Ri - ( Rj + T ) where 'i' core & 'j' LPcore (in a periodic image)
		// Evalulation
		Manager::support_h_matrix_real( lp, C.sigma, Rij, this->man_matrix4d_ws[0], this->man_matrix4d_ws[1] );

		factor = 0.5 * C.AtomList[i]->charge * lp->lp_charge;
		this->man_matrix4d_ws[1] = factor * this->man_matrix4d_ws[1];		// POSSIBLE MEMOIZATION ...

		// Calculate Energy Contribution by the given LP density in the periodic image and a point charge in the central sublattice
		partial_e = 0.;
		for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ partial_e += lp_cf[i] * lp_cf[j] * this->man_matrix4d_ws[1](i,j); }}			// PartialE ... Process Required
	}

	if( type_i == "lone" && type_j == "shel" )
	{
		LonePair* lp = static_cast<LonePair*>(C.AtomList[i]);

		// <1> W.R.T CorePart
		Rij = ( C.AtomList[j]->cart + TransVector ) - C.AtomList[i]->cart;	// (Rj+T) - Ri ... Ri(LPcore) ---> Rj+T(shell CorePart);
		// Evalulation
		Manager::support_h_matrix_real( lp, C.sigma, Rij, this->man_matrix4d_ws[0], this->man_matrix4d_ws[1] );

		factor = 0.5 * lp->lp_charge * C.AtomList[j]->charge;
		this->LPC_H_Real[i][j][0] += factor * this->man_matrix4d_ws[1];		// [0] for Core
		
		// <2> W.R.T ShelPart
		Rij = ( static_cast<Shell*>(C.AtomList[j])->shel_cart + TransVector ) - C.AtomList[i]->cart;	// (Rj+T) - Ri ... Ri(LPcore) ---> Rj+T(ShellPart);
		// Evalulation
		Manager::support_h_matrix_real( lp, C.sigma, Rij, this->man_matrix4d_ws[0], this->man_matrix4d_ws[1] );
		
		factor = 0.5 * lp->lp_charge * static_cast<Shell*>(C.AtomList[j])->shel_charge;
		this->LPC_H_Real[i][j][1] += factor * this->man_matrix4d_ws[1];		// [1] for Shell
	}

	if( type_i == "shel" && type_j == "lone" )
	{
		// Get 'j' lone pair cation & its eigenvectors of the ground-state
		LonePair* lp = static_cast<LonePair*>(C.AtomList[j]);
		lp_cf[0] = lp->lp_eigensolver.eigenvectors()(0,lp->lp_gs_index).real();	// s
		lp_cf[1] = lp->lp_eigensolver.eigenvectors()(1,lp->lp_gs_index).real();	// px
		lp_cf[2] = lp->lp_eigensolver.eigenvectors()(2,lp->lp_gs_index).real();	// py
		lp_cf[3] = lp->lp_eigensolver.eigenvectors()(3,lp->lp_gs_index).real();	// pz

		// <1> Core W.R.T LP Density 
		Rij = C.AtomList[i]->cart - ( C.AtomList[j]->cart + TransVector );	// Ri - ( Rj + T ) where 'i' core & 'j' LPcore (in a periodic image)
		// Evalulation
		Manager::support_h_matrix_real( lp, C.sigma, Rij, this->man_matrix4d_ws[0], this->man_matrix4d_ws[1] );

		factor = 0.5 * C.AtomList[i]->charge * lp->lp_charge;
		this->man_matrix4d_ws[1] = factor * this->man_matrix4d_ws[1];

		// Calculation Energy Contribution by the given LP density in the periodic image and the core in the central sublattice
		partial_e = 0.;
		for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ partial_e += lp_cf[i] * lp_cf[j] * this->man_matrix4d_ws[1](i,j); }}			// PartialE ... Process Required
		

		// <2> Shel W.R.T LP Density
		Rij = static_cast<Shell*>(C.AtomList[i])->shel_cart - ( C.AtomList[j]->cart + TransVector );
		// Evalulation
		Manager::support_h_matrix_real( lp, C.sigma, Rij, this->man_matrix4d_ws[0], this->man_matrix4d_ws[1] );

		factor = 0.5 * static_cast<Shell*>(C.AtomList[i])->shel_charge * lp->lp_charge;
		this->man_matrix4d_ws[1] = factor * this->man_matrix4d_ws[1];

		// Calculation Energy Contribution by the given LP density in the periodic image and the shel in the central sublattice
		partial_e = 0.;
		for(int i=0;i<4;i++){ for(int j=0;j<4;j++){ partial_e += lp_cf[i] * lp_cf[j] * this->man_matrix4d_ws[1](i,j); }}			// PartialE ... Process Required
	}

	if( type_i == "lone" && type_j == "lone" )
	{
		LonePair* lpi = static_cast<LonePair*>(C.AtomList[i]);
		LonePair* lpj = static_cast<LonePair*>(C.AtomList[j]);
	
		// <1> LP(i) Density ('central sublattice') vs LP(j) Core ('periodic image')
		Rij = ( C.AtomList[j]->cart + TransVector ) - C.AtomList[i]->cart;	// LP(i)('in the central sublattice')  -> LP(j)('in the periodic image') core
		// Evaluation
		Manager::support_h_matrix_real( lpi, C.sigma, Rij, this->man_matrix4d_ws[0], this->man_matrix4d_ws[1] );

		factor = 0.5 * lpi->lp_charge * C.AtomList[j]->charge;			// LP(i) charge * LP(j) core charge
		this->LPLP_H_Real[i][j] += factor * this->man_matrix4d_ws[1];

		// <2> LP(i) Core ('central sublattice') vs LP(j) Density ('periodic image')
		Rij = C.AtomList[i]->cart - ( C.AtomList[j]->cart + TransVector );	// LP(j)('in the periodic image') -> LP(j)('in the central sublattice')
		// Get LP(j) Density eigenvectors 
		lp_cf[0] = lpj->lp_eigensolver.eigenvectors()(0,lpj->lp_gs_index).real();	// s
		lp_cf[1] = lpj->lp_eigensolver.eigenvectors()(1,lpj->lp_gs_index).real();	// px
		lp_cf[2] = lpj->lp_eigensolver.eigenvectors()(2,lpj->lp_gs_index).real();	// py
		lp_cf[3] = lpj->lp_eigensolver.eigenvectors()(3,lpj->lp_gs_index).real();	// pz
		// Evalulation
		Manager::support_h_matrix_real( lpj, C.sigma, Rij, this->man_matrix4d_ws[0], this->man_matrix4d_ws[1] );

		factor = 0.5 * C.AtomList[i]->charge * lpj->lp_charge;			// LP(i) core charge * LP(j) charge

		// Calculation Energy Contribution by the given LP density in the periodic image and the LP core in the central sublattice
		partial_e = 0.;
		for(int i=0;i<4;i++){ for(int j=0;j<4;j++) { partial_e += lp_cf[i] * lp_cf[j] * this->man_matrix4d_ws[1](i,j); }}			// PartialE ... Process Required

		// <3> LP(i) Density vs LP(j) Density


	}


	return;
}

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
