#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <omp.h>

int main()
{
	double val = 345226033152000;
	double val2= 345226033152000;
	double val3=1467210640896000;
	double val4=1467210640896000;
	
	val=34519618525593600;
	std::cout << val << std::endl;
	
	double et = -omp_get_wtime();
	et += omp_get_wtime();

	std::cout << std::setprecision(20) << et << std::endl;

	return 0;
}
