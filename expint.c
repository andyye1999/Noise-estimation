#include "expint.h"
#include <stdio.h>
#include <stdlib.h>
struct cheb_series_struct {
    double* c;   /* coefficients                */
    int order;    /* order of expansion          */
    double a;     /* lower interval point        */
    double b;     /* upper interval point        */
    int order_sp; /* effective single precision order */
};
typedef struct cheb_series_struct cheb_series;

static double AE11_data[39] = {
   0.121503239716065790,
  -0.065088778513550150,
   0.004897651357459670,
  -0.000649237843027216,
   0.000093840434587471,
   0.000000420236380882,
  -0.000008113374735904,
   0.000002804247688663,
   0.000000056487164441,
  -0.000000344809174450,
   0.000000058209273578,
   0.000000038711426349,
  -0.000000012453235014,
  -0.000000005118504888,
   0.000000002148771527,
   0.000000000868459898,
  -0.000000000343650105,
  -0.000000000179796603,
   0.000000000047442060,
   0.000000000040423282,
  -0.000000000003543928,
  -0.000000000008853444,
  -0.000000000000960151,
   0.000000000001692921,
   0.000000000000607990,
  -0.000000000000224338,
  -0.000000000000200327,
  -0.000000000000006246,
   0.000000000000045571,
   0.000000000000016383,
  -0.000000000000005561,
  -0.000000000000006074,
  -0.000000000000000862,
   0.000000000000001223,
   0.000000000000000716,
  -0.000000000000000024,
  -0.000000000000000201,
  -0.000000000000000082,
   0.000000000000000017
};
static cheb_series AE11_cs = {
  AE11_data,
  38,
  -1, 1,
  20
};

static double AE12_data[25] = {
   0.582417495134726740,
  -0.158348850905782750,
  -0.006764275590323141,
   0.005125843950185725,
   0.000435232492169391,
  -0.000143613366305483,
  -0.000041801320556301,
  -0.000002713395758640,
   0.000001151381913647,
   0.000000420650022012,
   0.000000066581901391,
   0.000000000662143777,
  -0.000000002844104870,
  -0.000000000940724197,
  -0.000000000177476602,
  -0.000000000015830222,
   0.000000000002905732,
   0.000000000001769356,
   0.000000000000492735,
   0.000000000000093709,
   0.000000000000010707,
  -0.000000000000000537,
  -0.000000000000000716,
  -0.000000000000000244,
  -0.000000000000000058
};
static cheb_series AE12_cs = {
  AE12_data,
  24,
  -1, 1,
  15
};

static double E11_data[19] = {
  -16.11346165557149402600,
    7.79407277874268027690,
   -1.95540581886314195070,
    0.37337293866277945612,
   -0.05692503191092901938,
    0.00721107776966009185,
   -0.00078104901449841593,
    0.00007388093356262168,
   -0.00000620286187580820,
    0.00000046816002303176,
   -0.00000003209288853329,
    0.00000000201519974874,
   -0.00000000011673686816,
    0.00000000000627627066,
   -0.00000000000031481541,
    0.00000000000001479904,
   -0.00000000000000065457,
    0.00000000000000002733,
   -0.00000000000000000108
};
static cheb_series E11_cs = {
  E11_data,
  18,
  -1, 1,
  13
};

static double E12_data[16] = {
  -0.03739021479220279500,
   0.04272398606220957700,
  -0.13031820798497005440,
   0.01441912402469889073,
  -0.00134617078051068022,
   0.00010731029253063780,
  -0.00000742999951611943,
   0.00000045377325690753,
  -0.00000002476417211390,
   0.00000000122076581374,
  -0.00000000005485141480,
   0.00000000000226362142,
  -0.00000000000008635897,
   0.00000000000000306291,
  -0.00000000000000010148,
   0.00000000000000000315
};
static cheb_series E12_cs = {
  E12_data,
  15,
  -1, 1,
  10
};

static double AE13_data[25] = {
  -0.605773246640603460,
  -0.112535243483660900,
   0.013432266247902779,
  -0.001926845187381145,
   0.000309118337720603,
  -0.000053564132129618,
   0.000009827812880247,
  -0.000001885368984916,
   0.000000374943193568,
  -0.000000076823455870,
   0.000000016143270567,
  -0.000000003466802211,
   0.000000000758754209,
  -0.000000000168864333,
   0.000000000038145706,
  -0.000000000008733026,
   0.000000000002023672,
  -0.000000000000474132,
   0.000000000000112211,
  -0.000000000000026804,
   0.000000000000006457,
  -0.000000000000001568,
   0.000000000000000383,
  -0.000000000000000094,
   0.000000000000000023
};
static cheb_series AE13_cs = {
  AE13_data,
  24,
  -1, 1,
  15
};

static double AE14_data[26] = {
  -0.18929180007530170,
  -0.08648117855259871,
   0.00722410154374659,
  -0.00080975594575573,
   0.00010999134432661,
  -0.00001717332998937,
   0.00000298562751447,
  -0.00000056596491457,
   0.00000011526808397,
  -0.00000002495030440,
   0.00000000569232420,
  -0.00000000135995766,
   0.00000000033846628,
  -0.00000000008737853,
   0.00000000002331588,
  -0.00000000000641148,
   0.00000000000181224,
  -0.00000000000052538,
   0.00000000000015592,
  -0.00000000000004729,
   0.00000000000001463,
  -0.00000000000000461,
   0.00000000000000148,
  -0.00000000000000048,
   0.00000000000000016,
  -0.00000000000000005
};
static cheb_series AE14_cs = {
  AE14_data,
  25,
  -1, 1,
  13
};


static inline double cheb_eval(const cheb_series* cs,
    const double x)
{
    int j;
    double d = 0.0;
    double dd = 0.0;

    double y = (2.0 * x - cs->a - cs->b) / (cs->b - cs->a);
    double y2 = 2.0 * y;

    for (j = cs->order; j >= 1; j--)
    {
        double temp = d;
        d = y2 * d - dd + cs->c[j];
        dd = temp;
    }

    return y * d - dd + 0.5 * cs->c[0];

}

double expint_E1(double x, int scale) 
{
#ifdef IEEE_754
    if (ISNAN(x))
        return x;
#endif

    const double xmaxt = -LOG_DBL_MIN;       /* XMAXT = -LOG(DBL_MIN) */
    const double xmax = xmaxt - log(xmaxt); /* XMAX = XMAXT - LOG(XMAXT) */
    
    if (x < -xmax && !scale)
    {
        return;
    }

    else if (x <= -10.0)
    {
        const double s = 1.0 / x * (scale ? 1.0 : exp(-x));
        const double cheb = cheb_eval(&AE11_cs, 20.0 / x + 1.0);
        return s * (1.0 + cheb);
    }
    else if (x <= -4.0)
    {
        const double s = 1.0 / x * (scale ? 1.0 : exp(-x));
        const double cheb = cheb_eval(&AE12_cs, (40.0 / x + 7.0) / 3.0);
        return s * (1.0 + cheb);
    }
    else if (x <= -1.0)
    {
        const double s = (scale ? exp(x) : 1.0);
        const double ln_term = -log(fabs(x));
        const double cheb = cheb_eval(&E11_cs, (2.0 * x + 5.0) / 3.0);
        return s * (ln_term + cheb);
    }
    else if (x == 0.0)
    {
        return;
    }
    else if (x <= 1.0)
    {
        const double s = (scale ? exp(x) : 1.0);
        const double ln_term = -log(fabs(x));
        const double cheb = cheb_eval(&E12_cs, x);
        return s * (ln_term - 0.6875 + x + cheb);
    }
    else if (x <= 4.0)
    {
        const double s = 1.0 / x * (scale ? 1.0 : exp(-x));
        const double cheb = cheb_eval(&AE13_cs, (8.0 / x - 5.0) / 3.0);
        return s * (1.0 + cheb);
    }
    else if (x <= xmax || scale)
    {
        const double s = 1.0 / x * (scale ? 1.0 : exp(-x));
        const double cheb = cheb_eval(&AE14_cs, 8.0 / x - 1.0);
        double res = s * (1.0 + cheb);
        if (res == 0.0)
        {
            
            return 0.0;
        }
        else
            return res;
    }
    else {
        return 0;
    }
}

double expint_E2(double x, int scale) 
{
#ifdef IEEE_754
    if (ISNAN(x))
        return x;
#endif

    const double xmaxt = -LOG_DBL_MIN;
    const double xmax = xmaxt - log(xmaxt);
    if (x < -xmax && !scale)
    {
        return;
    }
    else if (x == 0.0)
    {
        return 1.0;
    }

    else if (x < 100.0)
    {
        const double ex = (scale ? 1.0 : exp(-x));
        return ex - x * expint_E1(x, scale);
    }
    else if (x < xmax || scale)
    {
        const double s = (scale ? 1.0 : exp(-x));
        const double c1 = -2.0;
        const double c2 = 6.0;
        const double c3 = -24.0;
        const double c4 = 120.0;
        const double c5 = -720.0;
        const double c6 = 5040.0;
        const double c7 = -40320.0;
        const double c8 = 362880.0;
        const double c9 = -3628800.0;
        const double c10 = 39916800.0;
        const double c11 = -479001600.0;
        const double c12 = 6227020800.0;
        const double c13 = -87178291200.0;
        const double y = 1.0 / x;
        const double sum6 = c6 + y * (c7 + y * (c8 + y * (c9 + y * (c10 + y * (c11 + y * (c12 + y * c13))))));
        const double sum = y * (c1 + y * (c2 + y * (c3 + y * (c4 + y * (c5 + y * sum6)))));
        double res = s * (1.0 + sum) / x;
        if (res == 0.0)
        {
            return 0.0;
        }
        else
            return res;
    }
    else {
        return 0.0;
    }
}

double gammainc(double s, double x, int p)
{
    int k = 0;
    double sum = 0;
    for (k = 0; k < p; k++) {
        sum += pow(x, k) / tgamma(s + k + 1);
    }
    return pow(x, s) * tgamma(s) * exp(-x) * sum;
}

void conv(float* Ptr_Src1, float* Ptr_Src2, int Src1Lenth, int Src2Lenth, float* Ptr_Target)
{
    float temp = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int ALL_Length;
    ALL_Length = Src1Lenth + Src2Lenth - 1;
    ///================卷积核心算法部分==============
    for (i = 0; i < ALL_Length; i++)
    {
        Ptr_Target[i] = 0;
    }
    for (i = 0; i < ALL_Length; i++)
    {
        for (j = max(0, i + 1 - Src2Lenth); j <= min(i, Src1Lenth - 1); j++)
        {
            Ptr_Target[i] += Ptr_Src1[j] * Ptr_Src2[i - j];
        }
        //printf("Ptr_Target[%d] =%f\n",i,Ptr_Target[i]);
        //getchar();
    }
}