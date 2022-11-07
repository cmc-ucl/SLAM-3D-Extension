#include <cstdio>
#include <cmath>
#include <gsl/gsl_integration.h>

double f (double x, void * params) {
    double alpha = *(double *) params;
    double f = alpha / (x * x + 1);
    return f;
}

    int
main (void)
{
    gsl_integration_workspace * w 
        = gsl_integration_workspace_alloc (1000);

    double result, error;
    double alpha = 1.0;


    gsl_function F;
    F.function = &f;
    F.params = &alpha;

    gsl_integration_qag (&F,
                         0.0, 1000.0,
                         0.0, 1e-7, 1000,
                         GSL_INTEG_GAUSS15,
                         w,
                         &result, &error);

    printf ("result          = % .18f\n", result);
    printf ("estimated error = % .18f\n", error);

    gsl_integration_workspace_free (w);

    return 0;
} 
