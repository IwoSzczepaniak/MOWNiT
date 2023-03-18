#include <stdio.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>

int main()
{
    gsl_ieee_env_setup();
    float number = 1e-34;

    while (number > 0)
    {
        number /= (2.0);
        gsl_ieee_printf_float(&number);
        printf("\n");
    }
    return 0;
}