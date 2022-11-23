#include <iostream>
#include <cmath>
#include <gsl/gsl_integration.h>

#include <omp.h>

double f (double x, void * params)
{
	double alpha = *(double *) params;
	double f = log(alpha*x) / sqrt(x);
	return f;
}

struct g_params{ double a; double b; double c;};


double g(double x, void* params)
{
	struct g_params* p = (struct g_params*)params;

	double a1 = p->a;
	double b1 = p->b;
	double sig= p->c;

	double g  = (a1*x*x + b1*x)*erf(x/sig);

	return g;
}

int main (void)
{
	double et;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

	double result, error;
	double expected = -4.0;
	double alpha = 1.0;

	gsl_function F;
	F.function = &f;
	F.params = &alpha;

	std::cout << "1> Range 0 - 1 -----------" << std::endl;

	et = -omp_get_wtime();
	gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000, w, &result, &error);
	et += omp_get_wtime();

	printf ("result          = % .18f\n", result);
	printf ("exact result    = % .18f\n", expected);
	printf ("estimated error = % .18f\n", error);
	printf ("actual error    = % .18f\n", result - expected);
	printf ("intervals       = %zu\n", w->size);
	printf ("et              = % .18f\n",et);


	std::cout << "2> Range 0.5 - 4 -----------" << std::endl;

	struct g_params p = {-2.4,0.9,1.123};

	gsl_function G;
	G.function = &g;
	G.params = &p;


	expected = -43.644530795660636;
	et = -omp_get_wtime();
	gsl_integration_qag(&G,0.0,4.0   ,0.0, 1e-7, 1000,GSL_INTEG_GAUSS15,w,&result, &error);
	et += omp_get_wtime();
	//gsl_integration_qag(&F,0.0,1000.0,0.0, 1e-7, 1000,GSL_INTEG_GAUSS15,w,&result, &error);

	printf ("result          = % .18f\n", result);
	printf ("exact result    = % .18f\n", expected);
	printf ("estimated error = % .18f\n", error);
	printf ("actual error    = % .18f\n", result - expected);
	printf ("intervals       = %zu\n", w->size);
	printf ("et              = % .18f\n",et);



	gsl_integration_workspace_free (w);

	return 0;
}
