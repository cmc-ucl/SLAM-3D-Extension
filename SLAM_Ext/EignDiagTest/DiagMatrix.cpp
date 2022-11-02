#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include <vector>

#include <fstream>
#include <sstream>

int main()
{
using std::cout, std::endl;


	cout << "# Test 1 .............................................." << endl;
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

	cout << "# Test 2 .............................................." << endl;
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

	cout << "Mid Evec" << endl;
	printf("%12.6lf\t%12.6lf\t%12.6lf\n",es.eigenvectors()(1,0).real(),es.eigenvectors()(1,1).real(),es.eigenvectors()(1,2).real());

	cout << "Get Min index" << endl;
	std::vector<double> evals;
	int min_idx;
	for(int i=0;i<3;i++)
	{	evals.push_back(es.eigenvalues()(i).real());
	}
	min_idx = std::min_element(evals.begin(),evals.end()) - evals.begin();
	cout << min_idx << endl;


	cout << endl;
	cout << endl;
	cout << endl;
	cout << endl;
	cout << "Testing Hessian" << endl;

	Eigen::MatrixXd H(12,12);

	std::ifstream fp;
	std::stringstream ss;
	fp.open("hess.txt");
	double tmp;

	int row = 0;
	for(std::string str; std::getline(fp,str); )
	{	
		ss.clear();
		ss.str("");
		ss << str;

		for(int j=0;j<12;j++)
		{
			ss >> tmp;
			printf("%12.6lf\t",tmp);
			H(row,j) = static_cast<double>(tmp);	
		}
		cout << endl;
		row++;
	}

	cout << "Print" << endl;
	for(int i=0;i<12;i++)
	{	for(int j=0;j<12;j++)
		{	printf("%18.12lf\t",H(i,j));
		}
		cout << endl;
	}

	cout << H(3,11) << endl;
	Eigen::EigenSolver<Eigen::MatrixXd> es2;		// WatchOut Template Type!	'< Eigen::MatrixXd >'
	cout << "\n\nDIAGONALISATION !!!" << endl;
/*
	//Eigen::MatrixXd t(4,4);
	Eigen::Matrix4d t;

	t << 1, 2, 0, 5,
	     2,-4, 3, 1,
	     0, 3,-1, 5,
	     5, 1, 5,10;
	cout << t << endl;
*/
	es2.compute(H,true);	// true flag : compute_eigenvectors=true
	cout << "Eval : \n";
	cout << es2.eigenvalues() << endl;
	cout << "EigenVectors : \n";
	cout << es2.eigenvectors() << endl;
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
