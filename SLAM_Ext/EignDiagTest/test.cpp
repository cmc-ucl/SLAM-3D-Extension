#include<iostream>
#include<Eigen/Dense>
/*
const Eigen::Matrix4d& foo( const double var, Eigen::Matrix4d& mat )
{
	mat.setZero();
	mat(0,0) = 1;
	return mat;
}
*/
void foo( const double var, Eigen::Matrix4d& mat )
{
	mat.setZero();
	mat(0,0) = var;
	
	mat(1,2) = mat(2,1) = var*10;
	
}

int main()
{
	using std::cout, std::endl;

	Eigen::Matrix4d m;
	Eigen::Matrix4d n;


	foo(10,m);

	cout << m << endl;

	cout << endl;
	cout << "Onsite operation" << endl;
	m = 1.23 * m;
	cout << m << endl;

	return 0;
}
