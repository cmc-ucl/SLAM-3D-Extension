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

void foo2( int (&arr1)[2], int (&arr2)[3] )
{
	std::cout << "1st Arr" << std::endl;
	for(auto x : arr1 )
	{
		std::cout << x << std::endl;
	}

	std::cout << "2nd Arr" << std::endl;
	for(auto x : arr2 )
	{
		std::cout << x << std::endl;
	}
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

	cout << "zero test" << endl;

	int h;

	h = 0;

	double val = h*0.213124112;

	if( val == 0 )
	{	cout << "yes " << endl;
	}	

	cout << "\n\n\n\n Test Passing ArrReference" << endl;
	int arr1[] = {4,5};
	int arr2[] = {1,2,3};	// SIZE MUST IN MATCH!!!

	foo2( arr1, arr2 );


	return 0;
}
