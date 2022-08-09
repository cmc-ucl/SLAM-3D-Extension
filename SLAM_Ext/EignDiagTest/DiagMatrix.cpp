#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include <vector>

int main()
{
using std::cout, std::endl;

	Eigen::EigenSolver<Eigen::Matrix3d> es;
	Eigen::Matrix3d m;

	m << 1.1,0.76,3.25,
	     0.76,4.6,5,
	     3.25,5,-6.1;

	cout << "1>Target Matrix: \n";
	cout << m << endl;
	es.compute(m,true);	// true flag : compute_eigenvectors=true
	cout << "Eval : \n";
	cout << es.eigenvalues() << endl;
	cout << "EigenVectors : \n";
	cout << es.eigenvectors() << endl;

	m << 1.1,-0.6,3.25,
	     -0.6,-4.6,5,
	     3.25,5,-6.1;

	cout << endl;
	cout << "2> Target Matrix : \n";
	cout << m << endl;
	es.compute(m,true);	// true flag : compute_eigenvectors=true
	cout << "Eval : \n";
	cout << es.eigenvalues() << endl;
	cout << "EigenVectors : \n";
	cout << es.eigenvectors() << endl;

	cout << "Handling\n";
	cout << es.eigenvalues()[0] << endl;		// decltype(es.eignevalues()[0]) == std::complex<double>
	cout << es.eigenvalues()[0].real() << endl;	// member function 'real()' of std::complex<T>

	cout << "Eigenvalues  are : " << endl;
	printf("%12.6lf%12.6lf%12.6lf\n",es.eigenvalues()[0].real(),es.eigenvalues()[1].real(),es.eigenvalues()[2].real());
	cout << "Eigenvectors are : " << endl; 

	for(int i=0;i<3;i++)
	{
		printf("%12.6lf%12.6lf%12.6lf\n",es.eigenvectors()(i,0).real(),es.eigenvectors()(i,1).real(),es.eigenvectors()(i,2).real());
	}

	cout << "Get Min index" << endl;
	std::vector<double> evals;
	int min_idx;
	for(int i=0;i<3;i++)
	{	evals.push_back(es.eigenvalues()(i).real());
	}
	min_idx = std::min_element(evals.begin(),evals.end()) - evals.begin();
	cout << min_idx << endl;

/*
std::vector<int> v = {5, 2, 8, 10, 9}; 
int maxElementIndex = std::max_element(v.begin(),v.end()) - v.begin();
int maxElement = *std::max_element(v.begin(), v.end());

int minElementIndex = std::min_element(v.begin(),v.end()) - v.begin();
int minElement = *std::min_element(v.begin(), v.end());

std::cout << "maxElementIndex:" << maxElementIndex << ", maxElement:" << maxElement << '\n';
std::cout << "minElementIndex:" << minElementIndex << ", minElement:" << minElement << '\n';
*/
	return 0;
}
