#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

double foo( const double x )
{
	return 3.*x*x + 2.*x + 1.;
}

int main(int argc, char* argv[])
{
	const int thc = strtol(argv[1],NULL,10);
	double a = 0, b = 3.;
	int div = 100000000;
	double dx = (b-a)/div;

	double partial_sum;
	double total_sum = 0.;

	double x1,x2;
	double h1,h2;
	double area;

	double elapsed_t;


	elapsed_t = -omp_get_wtime();
	#pragma omp parallel private(x1,x2,h1,h2) num_threads(thc) 	// Important >>> the private variables !!! could be in the race condition
	{								// 'a' , 'dx' are shared 
		#pragma omp for reduction(+:total_sum)
		for(int i=0;i<div-1;i++)
		{	
			x1 = a + dx*i;	// a, dx (constants) i (private) >>>  x1 should be private
			x2 = x1 + dx;	// x1 is private >>> x2 should be private

			h1 = foo(x1);	// x1, x2 privates >>>> h1, h2 privates
			h2 = foo(x2);

			//# pragma omp critical		// without reduction cluase ?
			total_sum += (h2+h1)*dx/2.;
		}
	}
	elapsed_t += omp_get_wtime();

	printf("Answer %lf\n",total_sum);	// should be ~ 39.000 // simple exe with time ... threading is much slower ...
	printf("threads in use : %d / elapsed time : %lf\n",thc,elapsed_t);


	return 0;
}
