#include <math.h>


#define LOG_DBL_MIN   (-7.0839641853226408e+02)
#define LOG_DBL_MAX    7.0978271289338397e+02
#define EULER_CNST     0.57721566490153286060651209008
#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))


double expint_E1(double x, int scale);
double expint_E2(double x, int scale);
void conv(float* Ptr_Src1, float* Ptr_Src2, int Src1Lenth, int Src2Lenth, float* Ptr_Target);